# -*- coding: utf-8 -*-
"""
    femagtools.asm
    ~~~~~~~~~~~~~~

    Reading ASM files

Slots Qr:

  Qr <= 1.25 Qs
  Qr != Qs +/- 2p
  Qr != Qs
  Qr != Qs +/- 1

"""
import pathlib
import logging
import logging.config
import numpy as np
import lmfit

logger = logging.getLogger('femagtools.asm')


def imcur(s, w1, u1, r1, ls1, lh, ls2, r2):
    """return currents i1r, i1i, i2r, i2i"""
    xs1 = w1*ls1
    xh = w1*lh
    x1 = xs1+xh
    xs2 = w1*ls2
    x2 = xs2+xh
    # solve voltage equations for steady state
    A = np.array((
        (r1, -x1, 0, -xh),
        (x1, r1, xh, 0),
        (0, -s*xh, r2, -s*x2),
        (s*xh, 0, s*x2, r2)))
    return np.linalg.solve(A, np.array((u1, 0, 0, 0)))


def torque(p, w1, u1, s, r1, ls1, lh, ls2, r2):
    """return torque"""
    if s == 0:
        return 0
    i2 = imcur(s, w1, u1, r1, ls1, lh, ls2, r2)[2:]
    return 3*p*r2/s/w1*(i2[0]**2 + i2[1]**2)


def fit_current(w1, u1, slip, r1, ls1, lh, i1, cosphi):
    def imcurr(s, ls2, r2):
        return np.array(
            [x + 1j*y
             for x, y in [imcur(sx, w1, u1, r1, ls1, lh, ls2, r2)[:2]
                          for sx in s]])
    model = lmfit.model.Model(imcurr)
    ls2_guess = 0.0
    r2_guess = np.mean(
        [u1/i1x*sx
         for i1x, sx in zip(i1, slip) if abs(i1x) > 1e-6])
    params = model.make_params(r2=r2_guess, ls2=ls2_guess)
    guess = lmfit.models.update_param_vals(params, model.prefix)
    i1c = np.array([x*pf - 1j*x*np.sqrt(1-pf**2)
                    for x, pf in zip(i1, cosphi)])
    r = model.fit(i1c, params=guess, s=slip, verbose=True)
    return r.params['ls2'].value, r.params['r2'].value


valmap = {
    'Stator windings voltage (RMS)[V]': 'u1',
    'Wdgs-connection: 0': 'wdgconn',
    'Nominal frequency': 'f1',
    'Stator phase winding resistamce': 'r1',
    'Stator phase winding resistance': 'r1',
    'Stator phase end-winding reactance': 'xs1',
    'Stator phase winding wires/slot side': 'num_wires',
    'Effect. air gap  length          [mm]': 'lfe',
    'Effect. rotor bar length (+end)  [mm]': 'rbarlen',
    'Number of Phases': 'num_phases',
    'Number of Pole pairs': 'p',
    'Number of Poles simulated': 'p_gen',
    'Number of parallel windings': 'num_par_wdgs',
    'Rotor Lamination                [kg]': 'lamweight',
    'Rotor Conductors                [kg]': 'conweight',
    'MC-File used in calculation': 'mcfile',
    'Losses[W/kg] in MC-File': 'felosscoeff',
    'Max. No. Iterations': 'maxiters',
    'Change of Perm. max %': 'permchg'}


def _read_sections(lines):
    """return list of ASM sections

    sections are either surrounded by lines starting with '[***'
    or by starting with any 'Input data' or 'Simulation Results'
    Args:
      param lines (list) lines of ASM file to read

    Returns:
      list of sections
    """

    section = []
    for line in lines:
        if ('[****' in line or
            'Input data' in line or
                'Simulation Results' in line):
            if section:
                # skip empty lines
                i = 0
                try:
                    while not section[i]:
                        i = i+1
                except IndexError:
                    i = i-1
                yield section[i:]
                section = []
        else:
            section.append(line.strip())
    yield section


def read_input_data(content):
    r = dict()
    for l in content:
        if '=' in l or ':' in l:
            d = '=' if '=' in l else ':'
            k, v = [s.strip() for s in l.split(d)[:2]]
            if k == 'Wdgs-connection: 0':
                c = int(float(l.split()[-1]))
                r[valmap[k]] = ['open', 'star', 'delta'][c]
            else:
                r[valmap[k]] = float(v)
        elif '\t' in l:
            k, v = [s.strip() for s in l.split('\t')[:2]]
            if valmap[k] == 'mcfile':
                r['mcfile'] = v
            else:
                r[valmap[k]] = float(v)

    for k in ('num_phases', 'p', 'p_gen',
              'num_par_wdgs', 'maxiters'):
        if k in r:
            r[k] = int(r[k])
    return r


def read_simulation_results(content):
    unit = 1e3
    resmap = {
        'Torque = P2/(s.omega) [Nm]': 'Tp2',
        'Rotor-Losses P2 [kW]': 'p2',
        'Rotor-Losses P2 [W]': 'p2',
        'Stator FE  Pfe1 [kW]': 'pfe1',
        'Stator FE  Pfe1 [W]': 'pfe1'}
    r = dict(s=[], T=[], u1=[], i1=[], p1=[], cosphi=[],
             f1=[], pfe1=[], p2=[], Tp2=[])
    for l in content:
        if l.startswith('S LIP'):
            if l.find('POWER[W]') > -1:
                unit = 1
                continue
        a = l.split()
        if len(a) == 7:
            for k, v in zip(r.keys(), a):
                r[k].append(float(v))
        elif a:
            a = l.split(':')[-1].split()
            if len(a) == 1:
                try:
                    k = resmap[l.split(':')[0]]
                    r[k].append(float(a[0]))
                except KeyError:
                    logger.warning('Key %s ignored', l.split(':')[0])
    if unit > 1:
        for k in ('p1', 'p2', 'pfe1'):
            r[k] = [x*unit for x in r[k]]
    r['s'] = [s/100 for s in r['s']]
    return r


def parident(w1, u1, i1, cosphi, s, r1, ls1, lh):
    """returns equivalent circuit parameters: r2, ls"""
    logger.info("w1 %s u1 %s i1 %s cosphi %s s %s r1 %s ls1 %s lh %s",
                w1, u1, i1, cosphi, s, r1, ls1, lh)
    ls2, r2 = fit_current(w1, u1, s,
                          r1, ls1, lh, i1, cosphi)

    xi = w1*(ls1+lh - lh**2/(ls2+lh))
    return dict(r2=r2, ls2=ls2, sk=r2/np.sqrt(xi**2 + r1**2))


def read(arg):
    """read asm file

        Args:
          filename or content (list of str) the text lines of the ASM file
    """
    r = {}
    if isinstance(arg, str):
        lines = pathlib.Path(arg).read_text().split('\n')
    elif isinstance(arg, pathlib.Path):
        lines = arg.read_text().split('\n')
    else:
        lines = arg
    for s in _read_sections(lines):
        if not s:
            continue
        title = s[0].split(':')[0].strip()

        if 'FEMAG Classic Version' in title:
            r['version'] = s[0].split(':')[-1].replace(' Version ', '').strip()
        elif 'Project File name' in title:
            r['project'] = s[1].strip()
        elif 'Number of Nodes' in title:
            pass
        elif 'File name' in title:
            r['filename'] = s[0].split(':')[-1].strip()
        elif 'Date' in title:
            d = s[0].split(':')[1].strip().split()
            dd, MM, yy = d[0].split('.')
            hh, mm = ''.join(d[1:-1]).split('.')
            r['date'] = '{}-{}-{}T{:02}:{:02}'.format(
                yy, MM, dd, int(hh), int(mm))
        elif 'Stator windings' in title:
            r.update(read_input_data(s))
        else:
            r.update(read_simulation_results(s))
    r['pcu'] = [3*x**2*r['r1'] for x in r['i1']]
    if r['pfe1']:
        r['pltotal'] = [sum(x) for x in zip(r['pcu'], r['pfe1'], r['p2'])]
    else:
        r['pltotal'] = [sum(x) for x in zip(r['pcu'], r['p2'])]

    w1 = 2*np.pi*r['f1'][0]
    r['ls1'] = r.pop('xs1')/w1
    u1 = r['u1'][0]
    i1 = r['i1']
    cosphi = r['cosphi']
    if r['wdgconn'] == 'star':
        u1 = u1/np.sqrt(3)
    if r['wdgconn'] == 'delta':
        i1 = [x/np.sqrt(3) for x in i1]

    r['lh'] = u1/i1[0]/w1 - r['ls1']
    if len(i1) > 2:
        r.update(parident(w1, u1, i1, cosphi,
                          r['s'], r['r1'], r['ls1'], r['lh']))
    return r


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    r = read(sys.argv[1])
    print(r)
