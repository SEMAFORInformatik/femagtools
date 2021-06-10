# -*- coding: utf-8 -*-
"""
    femagtools.asm
    ~~~~~~~~~~~~~~

    Reading ASM files



"""
import pathlib
import logging
import logging.config

logger = logging.getLogger('femagtools.asm')

valmap = {
    'Stator windings voltage (RMS)[V]': 'u1',
    'Wdgs-connection: 0': 'wdgconn',
    'Nominal frequency': 'f1',
    'Stator phase winding resistamce': 'r1',
    'Stator phase end-winding reactance': 'xs1',
    'Stator phase winding wires/slot side': 'num_wires',
    'Effect. air gap  length          [mm]': 'lfe',
    'Effect. rotor bar length (+end)  [mm]': 'rbarlen',
    'Number of Phases': 'num_phases',
    'Number of Pole pairs': 'p',
    'Number of Poles simulated': 'p_gen',
    'Number of parallel windings': 'num_par_wdgs',
    'MC-File used in calculation': 'mcfile',
    'Losses[W/kg] in MC-File': 'felosscoeff',
    'Max. No. Iterations': 'maxiters',
    'Change of Perm. max %': 'permchg'}


def _read_sections(f):
    """return list of ASM sections

    sections are either surrounded by lines starting with '[***'
    or by starting with any 'Input data' or 'Simulation Results'
    Args:
      param f (file) BCH file to be read

    Returns:
      list of sections
    """

    section = []
    for line in f:
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
            if 'MC-File' in k:
                r[valmap[k]] = v
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
        'Torque = P2/(s.omega) [Nm]': 'Tp',
        'Rotor-Losses P2 [kW]': 'p2',
        'Stator FE  Pfe1 [kW]': 'pfe1'}
    r = dict(s=[], T=[], un=[], i1=[], p1=[], cosphi=[],
             f1=[], pfe1=[], p2=[], Tp=[])
    for l in content:
        if l.startswith('S LIP'):
            if l.find('POWER[W]') > -1:
                unit = 1
                continue
        a = l.split()
        if len(a) == 7:
            for k, v in zip(('s', 'T', 'un', 'i1', 'p1', 'cosphi', 'f1'), a):
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
    return r


def read(arg):
    """read asm file

        Args:
          filename or content (list of str) the text lines of the ASM file
    """
    r = {}
    if isinstance(arg, str):
        lines = pathlib.Path(arg).read_text().split('\n')
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
    return r


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    r = read(argv[1])
    print(r)
