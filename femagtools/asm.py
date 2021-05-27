# -*- coding: utf-8 -*-
"""
    femagtools.asm
    ~~~~~~~~~~~~~~

    Reading ASM files



"""
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
    r = dict(s=[], T=[], un=[], i1=[], p1=[], cosphi=[], f1=[])
    for l in content:
        a = l.split()
        if len(a) == 7:
            for k, v in zip(('s', 'T', 'un', 'i1', 'p1', 'cosphi', 'f1'), a):
                r[k].append(float(v))
    return r


def read(content):
    """read asm file

        Args:
          content (str or list of str) the text lines of the ASM file
    """
    r = {}
    if isinstance(content, str):
        lines = content.split('\n')
    else:
        lines = content
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
    import pathlib
    import sys
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    content = pathlib.Path(sys.argv[1]).read_text()
    r = read(content)
    print(r)
