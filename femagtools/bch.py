# -*- coding: utf-8 -*-
"""
    femagtools.bch
    ~~~~~~~~~~~~~~

    Reading BCH/BATCH files



"""
import sys
import numpy as np
import re
import logging
import logging.config

logger = logging.getLogger('femagtools.bch')

alpha20 = 3.93e-3  # temperature coeff of copper

_indxOpPattern = re.compile(r'[a-zA-Z0-9_]+(\[[-0-9]+\])')
_statloss = re.compile(r'Fe-Losses\s*Stator')
_rotloss = re.compile(r'Fe-Losses\s*Rotor')


def splitindex(name):
    "returns listname and index if name contains index"
    m = _indxOpPattern.search(name)
    if m:
        return (name.split('[')[0],
                int(''.join(list(m.group(1))[1:-1])))
    return (name, None)


def floatnan(s):
    """converts string to float
    returns NaN on conversion error"""
    try:
        return float(s)
    except ValueError:
        return float('NaN')


def r1_20(r1, theta):
    return r1/(1+alpha20*(theta-20))


def get_si_factor(contentline):
    "extract the first pattern and return conversion factor for SI unit"
    pattern = re.compile(r"\[([A-Za-z/0-9]+)\]")
    search = pattern.search(contentline)
    if search:
        if search.group(1) in ('kW', 'kNm'):
            return 1e3
    return 1.0


def _readSections(f):
    """return list of bch sections

    sections are surrounded by lines starting with '[***'

    Args:
      param f (file) BCH file to be read

    Returns:
      list of sections
    """

    section = []
    for line in f:
        if line.startswith('[****'):
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


class Reader:
    """Reads a BCH/BATCH-File"""
    _numPattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')
    _numPatternNaN = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?|nan)\s*')

    def __init__(self):
        self._fft = None
        self.type = None
        self.filename = ''
        self.external_rotor = False
        self.project = ''
        self.date = ''
        self.version = ''
        # known types:
        # MULTIPLE CALCULATION OF FORCES AND FLUX
        # Fast cogging calculation OF FORCES AND FLUX
        # Fast LD-LQ-Identification
        # Fast Psid-Psiq-Identification
        # Fast PM-Synchronous-Motor Simulation
        # Characteristics of Permanent-Magnet-Motors
        # T(n) simulation from Ld-Lq-Psim values
        self.wdg = None
        self.wdgfactors = []
        self.torque = []
        self.torque_fft = []
        self.psidq = {}
        self.psidq_ldq = None
        self.machine = {}
        self.lossPar = {}
        self.flux = {}
        self.magnet = {}
        self.airgapInduction = {}
        self.flux_fft = {}
        self.linearForce = []
        self.linearForce_fft = []
        self.powerSituation = {}
        self.scData = {}
        self.dqPar = {}
        self.ldq = {}
        self.losses = []
        self.demag = []
        self.weights = []
        self.weight = {}
        self.characteristics = []
        self.areas = []
        self.current_angles = []
        self.dispatch = {
            'General Machine Data': Reader.__read_general_machine_data,
            'Weigths': Reader.__read_weights,
            'Number of Nodes': Reader.__read_nodes_and_mesh,
            'Windings Data': Reader.__read_dummy,
            'Winding-Factors': Reader.__read_winding_factors,
            'Machine Data': Reader.__read_machine_data,
            'CAD-Parameter Data': Reader.__read_dummy,
            'Torque-Force': Reader.__read_torque_force,
            'Airgap Induction Br': Reader.__read_airgapInduction,
            'Fourier Analysis': Reader.__read_fft,
            'Machine excitation': Reader.__read_dummy,
            'DQ-Parameter for open Winding Modell': Reader.__read_dq_parameter,
            'Magnet Data': Reader.__read_magnet_data,
            'Date': Reader.__read_date,
            'Inertia': Reader.__read_inertia,
            'Losses': Reader.__read_dummy,
            'Fe-Hysteresis- and Eddy current Losses[W] > 0.1 % Max':
            Reader.__read_hysteresis_eddy_current_losses,
            'Losses [W]': Reader.__read_losses,
            'Losses for speed [1/min]': Reader.__read_losses_tab,
            'Losses from PSID-Psiq-Identification for speed [1/min]':
            Reader.__read_losses_tab,
            'Magnet loss data': Reader.__read_magnet_loss_data,
            'Project File name': Reader.__read_project_filename,
            'File name': Reader.__read_filename,
            'Windings input data': Reader.__read_windings,
            'Control parameters for Loss calculation': Reader.__read_lossPar,
            'PSID-Psiq-Identification': Reader.__read_psidq,
            'Ld-Lq-Identifikation aus PSID-Psiq-Identification':
            Reader.__read_psidq_ldq,
            'Ld-Lq-Identification': Reader.__read_ldq,
            'Ld-Lq-Identification RMS-values': Reader.__read_ldq,
            'Machine Data Rotor': Reader.__read_dummy,
            'Current Angles defined from no-load test':
            Reader.__read_current_angles,
            'FEMAG Version': Reader.__read_version,
            'FEMAG Classic Version': Reader.__read_version,
            'Simulation Data': Reader.__read_simulation_data,
            'Area [mm**2]': Reader.__read_areas,
            'Basic Machine parameters': Reader.__read_dummy,
            'Winding': Reader.__read_dummy,
            'Calculation time [sec]': Reader.__read_calctime,
            'Results for Angle I-Up [Degree]': Reader.__read_dummy,
            'Demagnetisation': Reader.__read_demagnetization,
            'Transient short circuit': Reader.__read_short_circuit,
            'Peak winding currents': Reader.__read_peak_winding_currents,
            'Flux observed': Reader.__read_flux,
            'Linear Force': Reader.__read_linear_force,
            'Demagnetization Data': Reader.__read_demagnetization,
            'Power Situation (VZS)': Reader.__read_power_situation,
            '** Characteristics of Permanent-Magnet-Motors **':
            Reader.__read_characteristics}

    def getStep(self):
        """@returns displacement step of flux values"""
        if len(self.flux) > 0:
            return self.flux[0]['displ'][1]-self.flux[0]['displ'][0]
        return None

    def read(self, content):
        """read bch file

        Args:
          content (str or list of str) the text lines of the BCH file
        """
        if isinstance(content, str):
            lines = content.split('\n')
        else:
            lines = content
        for s in _readSections(lines):
            if not s:
                continue
            title = s[0].split(':')[0].strip()
            if title == 'Function':
                title = s[0].split(':')[1].strip()
            logger.debug("'%s': %d", title, len(s[1:]))
            # Check if we are finished with the fourier analysis part(s)
            if title != 'Fourier Analysis':
                self._fft = None
            else:
                title2 = s[0].split(':')[1].strip()
                for k in ['Airgap Induction Br']:
                    if k == title2[:len(k)]:
                        title = title2[:len(k)]
            if title in self.dispatch:
                self.dispatch[title](self, s)

        if len(self.weights) > 0:
            w = list(zip(*self.weights))
            self.weight['iron'] = sum(w[0])
            self.weight['conductor'] = sum(w[1])
            self.weight['magnet'] = sum(w[2])
            self.weight['total'] = sum([sum(l) for l in w])
        return self

    def __findNums(self, l):
        rec = self._numPattern.findall(l)
        if 3 * '*' in l:  # min 3 '*'
            li = re.sub(r'[\*]+', str(np.NaN), l)
            rec = self._numPatternNaN.findall(li)
            logger.debug("Numbers with **** In: %s", l)
            logger.debug("Out: %s", li)
            logger.debug("%d. NUMS: %s", len(rec), rec)
        return rec

    def __read_version(self, content):
        rec = content[0].split(':')
        self.version = rec[-1].replace(' Version ', '').strip()

    def __read_project_filename(self, content):
        self.project = content[1].strip()

    def __read_filename(self, content):
        self.filename = content[0].split(':')[1].strip()

    def __read_date(self, content):
        d = content[0].split(':')[1].strip().split()
        dd, MM, yy = d[0].split('.')
        hh, mm = ''.join(d[1:-1]).split('.')
        self.date = '{}-{}-{}T{:02}:{:02}'.format(yy, MM, dd, int(hh), int(mm))

    def __read_nodes_and_mesh(self, content):
        self.nodes, self.elements, self.quality = \
            [floatnan(r[0]) for r in [self.__findNums(l)
                                      for l in content[:3]]]
        for l in content[3:]:
            m = re.match(r'\*+([^\*]+)\*+', l)
            if m:
                self.type = m.group(1).strip()
                return

    def __read_magnet_data(self, content):
        """read magnet data section"""
        for l in content:
            for v in [["Remanence Br [T]", 'Br'],
                      ["Co. Force Hc", 'Hc'],
                      ["Rel. Permeability", 'muer'],
                      ["Magn.Temperatur[grad]", 'Tmag'],
                      ["Temp.coeficient Br", 'alpha'],
                      ["Limit demagnetization   [%]", 'demag_pc'],
                      ["Limit demagnetization Hx [kA/m]", 'demag_hx'],
                      ["Total magnet area", 'area']]:

                if l.startswith(v[0]):
                    rec = self.__findNums(l)
                    if len(rec) > 0:
                        self.magnet[v[1]] = floatnan(rec[-1])
                    break

    def __read_windings(self, content):
        "read winding section"
        index = False
        wdgs = []
        for line in content:
            rec = line.split()
            if len(rec) > 0 and not index:
                index = rec[0].startswith('Index')
            elif len(rec) == 5 and index:
                if rec[0].startswith('Number'):
                    index = False
                else:
                    wdgs.append((abs(int(rec[1])),
                                 1 if int(rec[1]) > 0 else -1,
                                 float(rec[2]), float(rec[3]),
                                 float(rec[4])))
            else:
                index = False
        phases = set(list(zip(*wdgs))[0])
        self.windings = dict([(k, dict([(j, [])
                                        for j in ('dir',
                                                  'N',
                                                  'R',
                                                  'PHI')]))
                              for k in phases])
        for w in wdgs:
            self.windings[w[0]]['dir'].append(w[1])
            self.windings[w[0]]['N'].append(w[2])
            self.windings[w[0]]['R'].append(w[3]/1e3)
            self.windings[w[0]]['PHI'].append(w[4])

    def __read_winding_factors(self, content):
        "read winding factors section"
        self.wdgfactors = []
        for line in content:
            if line.startswith('Winding-Key'):
                self.wdgfactors.append(dict(order=[], wfac=[],
                                            skewf=[], total=[]))
            else:
                rec = self.__findNums(line)
                if len(rec) == 4:
                    self.wdgfactors[-1]['order'].append(int(rec[0]))
                    self.wdgfactors[-1]['wfac'].append(float(rec[1]))
                    self.wdgfactors[-1]['skewf'].append(float(rec[2]))
                    self.wdgfactors[-1]['total'].append(float(rec[3]))

    def __read_dummy(self, content):
        return

    def __read_magnet_loss_data(self, content):
        for line in content:
            if line.startswith('El.Conductivity Magnet'):
                self.magnet['sigma_PM'] = float(line.split()[-1])

    def __read_calctime(self, content):
        try:
            self.calctime = float(content[1])
        except (IndexError, ValueError):
            pass

    def __read_simulation_data(self, content):
        for i, line in enumerate(content):
            if line.startswith('Number of Phases m'):
                self.machine['m'] = int(float((line.split()[-1])))
            elif (line.startswith('Radius air-gap center (torque)') or
                  line.startswith('Air-gap center 1:radius/position')):
                rec = [l.strip() for l in line.split()]
                fc_radius = 1e-3*float(rec[-1])
                self.machine['fc_radius'] = fc_radius
            elif line.startswith('Number of parallel Windings'):
                self.machine['num_par_wdgs'] = int(float(line.split()[-1]))
            elif line.startswith('leak_dist_wind'):
                self.leak_dist_wind = self.__read_leak_dist_wind(content[i+1:])
        return

    def __read_leak_dist_wind(self, content):
        leak_wind = dict()
        vmap = {
            'Number of segments': 'nseg',
            'Number of poles simulated': 'npolsim',
            'Radius air-gap center (torque) [mm]': 'fc_radius',
            'Radius of perimeter [mm]': 'perimrad',
            'Bending radius vertical [mm]': 'vbendrad',
            'End winding height [mm]': 'endheight',
            'Wire diameter [mm]': 'wiredia',
            'Armature length [mm]': 'armatureLength',
            'End winding length [mm]': 'wdgendlen',
            'External inductance L0e/winding [H]': 'L0e',
            'Lde/winding [H]': 'Lde',
            'Lqe/winding [H]': 'Lqe',
            'Internal inductance L0i/winding [H]': 'L0i',
            'Ldi/winding [H]': 'Ldi',
            'Lqi/winding [H]': 'Lqi'}

        for line in content:
            rec = [l.strip() for l in line.split()]
            try:
                k = vmap[' '.join(rec[:-1])]
                v = float(rec[-1])
                if rec[-2] == '[mm]':
                    leak_wind[k] = 1e-3*v
                elif rec[-2] == '[H]':
                    leak_wind[k] = v
                else:
                    leak_wind[k] = int(v)
            except (KeyError, ValueError):
                pass
        return leak_wind

    def __read_current_angles(self, content):
        self.current_angles = []
        for l in content:
            rec = self.__findNums(l)
            if len(rec) == 3:
                self.current_angles.append(floatnan(rec[-1]))
        return

    def __read_lossPar(self, content):
        self.lossPar = {
            'fo': [],
            'Bo': [],
            'ch': [],
            'cw': [],
            'hf': [],
            'ef': [],
            'ic': [],
            'gamfe': [],
            'thetaw': [],
            'fillfe': [],
            'culength': [],
            'mat': []}
        for l in content:
            for v in [["Base Frequency", 'fo'],
                      ["Base Induction", 'Bo'],
                      ["Hysteresis-Coefficient", 'ch'],
                      ["Eddycurrent-Coefficient", 'cw'],
                      ["Hysteresis-Frequency-Coefficient", 'hf'],
                      ["Eddycurrent-Frequency-Coefficient", 'ef'],
                      ["Induction-Coefficient", 'ic'],
                      ["Specific Weight Iron", 'gamfe'],
                      ["Conductor Temperature", 'thetaw'],
                      ["Fillfactor Iron", 'fillfe'],
                      ["Relative conductor length", 'culength'],
                      ["Material factor", 'mat']]:

                if l.find(v[0]) > -1:
                    rec = self.__findNums(l)
                    if len(rec) > 0:
                        if v[1] == 'culength':
                            self.lossPar[v[1]].append(floatnan(rec[-1])/100)
                        else:
                            self.lossPar[v[1]].append(floatnan(rec[-1]))
                    break

        return l

    def __read_demagnetization(self, content):
        keys = [('displ', 'current_1', 'current_2', 'current_3',
                 'H_max', 'H_av', 'area'),
                ('displ', 'current', 'beta',
                 'H_max', 'H_av', 'area')]

        demag = dict()
        for l in content:
            rec = self.__findNums(l)
            if len(rec) == 7 or len(rec) == 6:
                i = 0 if len(rec) == 7 else 1
                self.demag.append(
                    {k: floatnan(r) for r, k in zip(rec, keys[i])})

                if demag and 'segment' in demag:
                    self.demag[-1].update(demag)
            else:
                for v in [["Limit Hc value", "lim_hc"],
                          ["Max. Magnetization", "br_max"],
                          ["Demagnetisation", "segment"],
                          ["Min. Magnetization", "br_min"],
                          ["Area demagnetized", "area"]]:
                    if l.find(v[0]) > -1:
                        if len(rec) > 0:
                            if v[1] == 'segment':
                                demag[v[1]] = int(rec[-1])
                            else:
                                demag[v[1]] = floatnan(rec[-1])
                        break
        if demag:
            if not 'segment' in demag:
                self.demag.append(demag)

    def __read_short_circuit(self, content):
        "read short circuit section"
        if content[2].startswith('Time'):
            m = []
            for l in content:
                rec = l.split()
                if len(rec) == 5 and not rec[0].startswith('Time'):
                    m.append([floatnan(x) for x in rec])

            m = np.array(m).T
            self.scData['time'] = m[0].tolist()
            self.scData['ia'] = m[1].tolist()
            self.scData['ib'] = m[2].tolist()
            self.scData['ic'] = m[3].tolist()
            self.scData['torque'] = m[4].tolist()
            return

        l = content[-1]
        rec = l.split()
        self.scData['speed'] = floatnan(rec[0])
        self.scData['ikd'] = floatnan(rec[1])
        self.scData['tkd'] = floatnan(rec[2])
        self.scData['iks'] = floatnan(rec[3])
        self.scData['tks'] = floatnan(rec[4])

    def __read_peak_winding_currents(self, content):
        self.scData['peakWindingCurrents'] = [float(x)
                                              for x in re.findall(r'[-0-9.]+',
                                                                  ''.join(content))]

    def __read_general_machine_data(self, content):
        for l in content:
            if l.find('Armature Length [mm]:') > -1:
                self.armatureLength = floatnan(l.split()[-1])
            elif l.find('Magn. Fluss Psim RMS') > 0:
                self.machine['psim'] = floatnan(l.split()[-1])
            elif l.find('Number of Pole pairs') > -1:
                self.machine['p'] = int(l.split()[-1])
            elif l.find('Number of Poles simulated') > -1:
                self.machine['p_sim'] = int(l.split()[-1])
            elif l.find('Total Number of Slots') > -1:
                self.machine['Q'] = int(l.split()[-1])
            elif l.find('Number of Slot-Sides sim.') > -1:
                self.machine['qs_sim'] = int(l.split()[-1])
            elif l.find('POC-File used in calculation') > -1:
                self.machine['pocfile'] = l.split(
                    ':')[-1].strip().replace('\\', '\\\\')
            elif l.find('MC-File used in calculation') > -1:
                self.machine['mcfile'] = l.split(
                )[-1].strip().replace('\\', '\\\\')
            elif l.find('Rotation Fraction') > -1:
                self.machine['period_frac'] = int(l.split()[-1])

    def __read_characteristics(self, content):
        characteristics = {}
        for i, l in enumerate(content):
            if l.startswith('[[***'):
                break
            for v in [['Voltage (operat. limit)', 'u1nom'],
                      ['Current (operat. limit)', 'i1nom'],
                      ['Angle I VS Up (speed = 0)', 'beta0'],
                      ['Resistance stator winding', 'r1'],
                      ['Inductance(I,Angle I-Up)', 'Ldnom'],
                      ['Inductance(I,Angle I-Up)', 'Lqnom'],
                      ['Power (operating limit)',  'Pnom'],
                      ['Stator ewdg inductance',   'Le'],
                      ['Stator external inductance Lex', 'Lex'],
                      ['Magn. flux (RMS)', 'psimnom'],
                      ['Effect. armature length', 'lfe'],
                      ['Power cut-off speed NC', 'nc'],
                      ['Number of Pole pairs', 'p'],
                      ['Max. current (RMS)', 'i1max'],
                      ['Rel. number wdg turns(wdg.1)', 'relw'],
                      ['Number of Phases', 'm'],
                      ['Min. speed', 'nmin'],
                      ['Max. speed', 'nmax']]:

                if l.find(v[0]) > -1:
                    rec = self.__findNums(l)
                    if len(rec) > 0:
                        characteristics[v[1]] = floatnan(rec[-1])
                        break

        characteristics['ldq'] = {}
        m = []
        if content[i+1].startswith('Wdg'):
            for k, l in enumerate(content[i+3:]):
                if l.startswith('[[***'):
                    break
                rec = l.split('\t')
                if len(rec) == 6:
                    m.append([floatnan(x) for x in rec])
        else:
            k = -3
        if m:
            m = np.array(m).T
            ncols = len(set(m[1]))
            i1 = np.reshape(m[0], (-1, ncols)).T[0]
            nrows = len(i1)

            logger.info('characteristics ld-lq %d x %d', nrows, ncols)
            characteristics['ldq'] = {
                'beta': m[1][:ncols][::-1].tolist(),
                'i1': i1.tolist(),
                'ld': (characteristics['lfe']*np.reshape(
                    m[2], (nrows, ncols)).T[::-1]).tolist(),
                'lq': (characteristics['lfe']*np.reshape(
                    m[3], (nrows, ncols)).T[::-1]).tolist(),
                'psim': (characteristics['lfe']*np.reshape(
                    m[4], (nrows, ncols)).T[::-1]).tolist(),
                'torque': (characteristics['lfe']*np.reshape(
                    m[5], (nrows, ncols)).T[::-1]).tolist()}

        m = []
        columns = [['n', 'id', 'iq', 'torque', 'p2'],
                   ['beta', 'cos_phi', 'u1', 'um'],
                   ['lang', 'ud', 'uq', 'i1'],
                   ['lang', 'ld', 'lq', 'psim']]
        nsec = 0
        characteristics['speed_torque'] = {}
        for l in content[k+i+6:]:
            if l.startswith('[[***'):
                break
            if not l:
                continue
            if l.startswith('Speed') and m:
                if nsec == 0:
                    m = np.array(m).T
                else:
                    m = np.array(m).T[1:]
                for j, k in enumerate(columns[nsec]):
                    characteristics['speed_torque'][k] = m[j].tolist()
                m = []
                nsec += 1
            else:
                rec = self.__findNums(l)
                if len(rec) > 3:
                    m.append([floatnan(x) for x in rec])

        self.characteristics.append(characteristics)

    def __read_flux(self, content):
        "read and append flux section"

        f = {'displ': [], 'flux_k': [], 'voltage_dpsi': [],
             'voltage_four': [], 'current_k': [], 'voltage_ir': []}

        for l in content:
            rec = l.split()
            if l.startswith('Flux-Area'):
                areas = self.__findNums(l)
                if not areas:
                    continue

                self.wdg = areas[0] if len(areas) == 1 else '{}-{}'.format(
                    areas[0], areas[1])
                if self.wdg not in self.flux:
                    self.flux[self.wdg] = []
            elif len(rec) == 7:
                f['displ'].append(floatnan(rec[1].strip()))
                f['flux_k'].append(floatnan(rec[2].strip()))
                f['voltage_dpsi'].append(floatnan(rec[3].strip()))
                f['voltage_four'].append(floatnan(rec[4].strip()))
                f['current_k'].append(floatnan(rec[5].strip()))
                f['voltage_ir'].append(floatnan(rec[6].strip()))
            elif rec and rec[0].startswith('['):
                f['displunit'] = re.search(r"\[([^\]]*)\]", l).group(1).strip()

        self.flux[self.wdg].append(f)
        self._fft = Reader.__read_flux_fft

    def __read_linear_force(self, content):
        "read and append linear force section"
        cosys = 'xy'
        f = {'displ': [], 'magnet_1': [], 'force_x': [],
             'force_y': [], 'f_idpsi': []}

        for l in content:
            rec = self.__findNums(l)
            if len(rec) > 4:
                f['displ'].append(floatnan(rec[1].strip()))
                f['magnet_1'].append(floatnan(rec[2].strip()))
                f['force_x'].append(floatnan(rec[3].strip()))
                f['force_y'].append(floatnan(rec[4].strip()))
                # TODO f['f_idpsi'].append(floatnan(rec[5].strip()))
            elif l.split()[-1] == 'Force_Z':
                cosys = 'rz'

        if cosys == 'rz':
            f['force_r'] = f.pop('force_x')
            f['force_z'] = f.pop('force_y')

        if len(f['displ']) > 0:
            if cosys == 'xy':
                ripple = [max(f['force_x']) - min(f['force_x']),
                          max(f['force_y']) - min(f['force_y'])]
                f['ripple_x'] = ripple[0]
                f['ripple_y'] = ripple[1]
            else:
                ripple = [max(f['force_r']) - min(f['force_r']),
                          max(f['force_z']) - min(f['force_z'])]
                f['ripple_r'] = ripple[0]
                f['ripple_z'] = ripple[1]

            self.linearForce.append(f)

            self._fft = Reader.__read_linearForce_fft

    def __read_linearForce_fft(self, content):
        "read and append linear force fft section"
        if not self._fft:
            return
        linearForce_fft = dict(order=[], force=[], force_perc=[],
                               a=[], b=[])
        for l in content:
            rec = self.__findNums(l)
            if len(rec) > 2:
                linearForce_fft['order'].append(int(rec[0].strip()))
                linearForce_fft['force'].append(floatnan(rec[1].strip()))
                linearForce_fft['force_perc'].append(floatnan(rec[2].strip()))
                if len(rec) > 4:
                    linearForce_fft['a'].append(floatnan(rec[3].strip()))
                    linearForce_fft['b'].append(floatnan(rec[4].strip()))
                else:
                    linearForce_fft['a'].append(0.0)
                    linearForce_fft['b'].append(0.0)

        if linearForce_fft['order']:
            self.linearForce_fft.append(linearForce_fft)

    def __read_fft(self, content):
        if self._fft:
            self._fft(self, content)

    def __read_flux_fft(self, content):
        "read and append flux fft section"

        flux_fft = dict(order=[], flux=[], flux_perc=[],
                        voltage=[], voltage_perc=[], a=[], b=[])
        for l in content:
            rec = self.__findNums(l)
            if len(rec) > 4:
                flux_fft['order'].append(int(rec[0].strip()))
                flux_fft['flux'].append(floatnan(rec[1].strip()))
                flux_fft['flux_perc'].append(floatnan(rec[2].strip()))
                flux_fft['voltage'].append(floatnan(rec[3].strip()))
                flux_fft['voltage_perc'].append(floatnan(rec[4].strip()))
                if len(rec) > 3:
                    flux_fft['a'].append(floatnan(rec[3].strip())*1e-2)
                    flux_fft['b'].append(floatnan(rec[4].strip())*1e-2)
                else:
                    flux_fft['a'].append(0.0)
                    flux_fft['b'].append(0.0)
        if self.wdg not in self.flux_fft:
            self.flux_fft[self.wdg] = []
        self.flux_fft[self.wdg].append(flux_fft)

    def __read_airgapInduction(self, content):
        "read and append airgapInduction section"
        import scipy.integrate as si
        import math

        logger.debug('read airgapInduction')
        i1beta = False  # format is either i1/beta or id/iq
        if 'i1' in self.ldq and 'beta' in self.ldq:
            i1 = self.ldq['i1']
            beta = self.ldq['beta']
        elif 'id' in self.psidq and 'iq' in self.psidq:
            id = self.psidq['id']
            iq = self.psidq['iq']
        else:
            i1 = []
            beta = []
            id = []
            iq = []

        an = [[], [], [], []]
        bn = [[], [], [], []]
        Bm = []
        Ba = []

        for line in content[5:]:
            if line.startswith('[****'):
                break
            if line.startswith("Current"):
                i1beta = True
                continue
            if line.startswith("C_STEP"):
                return   # ignore this section
            try:
                rec = self.__findNums(line)
                if len(rec) == 10:
                    f = [float(s) for s in rec]
                    an[0].append(f[2])
                    bn[0].append(f[3])
                    an[1].append(f[4])
                    bn[1].append(f[5])
                    an[2].append(f[6])
                    bn[2].append(f[7])
                    an[3].append(f[8])
                    bn[3].append(f[9])

                    a = (an[0][-1], an[1][-1], an[2][-1], an[3][-1])
                    b = (bn[0][-1], bn[1][-1], bn[2][-1], bn[3][-1])

                    def B(x):
                        return sum(a[i] * np.cos((2 * i + 1) * x) +
                                   b[i] * np.sin((2 * i + 1) * x)
                                   for i in (0, 1, 2, 3) if not math.isnan(a[i]) and not math.isnan(b[i]))

                    def Bdc(x):
                        return abs(B(x))

                    Ba.append(si.quad(Bdc, 0, 2 * np.pi,
                                      limit=250)[0] / (2 * np.pi))
                    Bm.append(max([B(x) for x in np.linspace(
                        0, 2 * np.pi, 100)]))

            except Exception as e:
                logger.debug("Conversion error: {} :: {}".format(e, line))

        self.airgapInduction = dict()

        if i1beta:
            ncols = len(beta)
            if ncols:
                nrows = len(Ba)//ncols
            else:
                nrows = len(i1)
            self.airgapInduction['beta'] = beta
            self.airgapInduction['i1'] = i1[:nrows]
        else:
            ncols = len(iq)
            self.airgapInduction['iq'] = iq
            self.airgapInduction['id'] = id
            nrows = len(self.airgapInduction['id'])

        try:
            self.airgapInduction['an'] = [np.reshape(an[j][:nrows*ncols],
                                                     (nrows, ncols)).T.tolist()
                                          for j in (0, 1, 2, 3)]
            self.airgapInduction['bn'] = [np.reshape(bn[j][:nrows*ncols],
                                                     (nrows, ncols)).T.tolist()
                                          for j in (0, 1, 2, 3)]
            self.airgapInduction['Bm'] = np.reshape(Bm[:nrows*ncols],
                                                    (nrows, ncols)).T.tolist()
            self.airgapInduction['Ba'] = np.reshape(Ba[:nrows*ncols],
                                                    (nrows, ncols)).T.tolist()
            # check for nan:
            if len(self.airgapInduction['an'][0]) > 1 and \
               len(self.airgapInduction['an'][0][0]) != len(self.airgapInduction['an'][0][1]):
                self.airgapInduction['an'] = [self.airgapInduction['an'][i][1:]
                                              for i in range(3)]
                self.airgapInduction['bn'] = [self.airgapInduction['bn'][i][1:]
                                              for i in range(3)]
                self.airgapInduction['Ba'] = self.airgapInduction['Ba'][1:]
                self.airgapInduction['Bm'] = self.airgapInduction['Bm'][1:]

                if len(self.airgapInduction['an'][0]) > 1 and \
                   len(list(filter(lambda x: np.isnan(x),
                                   list(zip(*self.airgapInduction['an'][0]))[0]))) > 0:
                    self.airgapInduction['an'] = [self.airgapInduction['an'][i][1:]
                                                  for i in range(3)]
                    self.airgapInduction['bn'] = [self.airgapInduction['bn'][i][1:]
                                                  for i in range(3)]
                    self.airgapInduction['Ba'] = zip(*zip(*self.airgapInduction['Ba'])
                                                     [1:])
                    self.airgapInduction['Bm'] = zip(*zip(*self.airgapInduction['Bm'])
                                                     [1:])
        except ValueError:
            print(self.airgapInduction['i1'])

    def __read_torque_force(self, content):
        "read and append force/torque section"

        torque = {
            'angle': [],
            'current_1': [],
            'force_x': [],
            'force_y': [],
            't_idpsi': [],
            'torque': []}
        for l in content:
            rec = self.__findNums(l)
            if len(rec) == 7:
                torque['angle'].append(floatnan(rec[1].strip()))
                torque['current_1'].append(floatnan(rec[2].strip()))
                torque['force_x'].append(floatnan(rec[3].strip()))
                torque['force_y'].append(floatnan(rec[4].strip()))
                torque['t_idpsi'].append(floatnan(rec[5].strip()))
                torque['torque'].append(floatnan(rec[6].strip()))

        if len(torque['angle']) > 0:
            ripple = max(torque['torque']) - min(torque['torque'])
            torque['ripple'] = ripple
            self.torque.append(torque)
        self._fft = Reader.__read_torque_force_fft

    def __read_torque_force_fft(self, content):
        "read and append force/torque fft section"
        columns = content[3].split()
        if len(columns) > 1 and columns[1] == 'Torque':
            torque_fft = dict(order=[],
                              torque=[],
                              torque_perc=[],
                              a=[],
                              b=[])
            for l in content:
                rec = self.__findNums(l)
                if len(rec) > 2:
                    torque_fft['order'].append(int(rec[0].strip()))
                    torque_fft['torque'].append(floatnan(rec[1].strip()))
                    torque_fft['torque_perc'].append(floatnan(rec[2].strip()))
                    if len(rec) > 3:
                        torque_fft['a'].append(floatnan(rec[3].strip())*1e-2)
                        torque_fft['b'].append(floatnan(rec[4].strip())*1e-2)
                    else:
                        torque_fft['a'].append(0.0)
                        torque_fft['b'].append(0.0)
            if torque_fft['order']:
                self.torque_fft.append(torque_fft)

    def __read_power_situation(self, content):
        "read and append power situation section"
        ps = dict()
        beta = ''
        for l in content:
            rec = self.__findNums(l)
            if len(rec) == 7:
                w = rec[0]
                ps[w] = dict(
                    voltage=floatnan(rec[1].strip()),
                    current_1=floatnan(rec[2].strip()),
                    beta=floatnan(rec[3].strip()),
                    cosphi=floatnan(rec[4].strip()),
                    powerp=1e3*floatnan(rec[5].strip()),
                    powerq=1e3*floatnan(rec[6].strip()))
            elif l.startswith('Angle current'):
                beta = rec[0]
        if ps:
            self.powerSituation = dict(windings=ps)
            if beta:
                self.powerSituation['beta'] = float(beta)

    def __removeTrailingZero(self, idList):
        '''if id list is inhomogeneous, remove a trailing '0'
        e.g. : idList[-450, -350, -250, -150, -50, 0]
               idList[-500, -400, -300, -200, -100, 0, 0]
        '''
        if idList[-1] == 0 and len(idList) > 2 and \
           int(idList[-1] - idList[0]) / (len(idList)-1) != \
           int(idList[-2] - idList[0]) / (len(idList)-2):
            idList = idList[:-1]
        return idList

    def __read_psidq(self, content):
        "read psid-psiq section"
        for i, l in enumerate(content):
            if l.find('[A') > -1:
                break
        m = []
        for l in content[i+2:]:
            rec = l.split('\t')
            if len(rec) == 7:
                m.append([floatnan(x) for x in rec])

        m = np.array(m).T
        ncols = np.argmax(np.abs(m[1][1:]-m[1][:-1]))+1
        if ncols == 1 and len(m[1]) > 1 and m[1][0] != m[1][1]:  # simple correction
            ncols = 2
        iq = np.reshape(m[1], (-1, ncols))[0]
        if ncols > 1 and (iq[ncols-1] < iq[ncols-2] or
                          len(m[0]) % ncols != 0):
            ncols = ncols-1

        id = np.reshape(m[0], (-1, ncols)).T[0]
        id = self.__removeTrailingZero(id)
        nrows = len(id)
        if nrows > 1 and id[nrows-1] < id[nrows-2]:
            nrows = nrows-1
        logger.info('psid-psiq %d x %d', nrows, ncols)
        mlen = nrows*ncols
        self.psidq = {k: (self.armatureLength*np.reshape(
            v[:mlen], (nrows, ncols))).T.tolist()
            for k, v in zip(('psid', 'psiq', 'torque_fe', 'torque'),
                            m[3:])}
        self.psidq['iq'] = iq[:ncols].tolist()
        self.psidq['id'] = id[:nrows].tolist()

    def __read_psidq_ldq(self, content):
        "read ldq from psid-psiq section"
        for i, l in enumerate(content):
            if l.find('[A') > -1:
                break
        m = []
        for l in content[i+2:]:
            rec = l.split('\t')
            if len(rec) == 7:
                m.append([floatnan(x) for x in rec])

        m = np.array(m).T
        ncols = np.argmax(np.abs(m[1][1:]-m[1][:-1]))+1
        if ncols == 1 and len(m[1]) > 1 and m[1][0] != m[1][1]:  # simple correction
            ncols = 2
        iq = np.linspace(np.min(m[1]), np.max(m[1]), ncols)
        if ncols > 1 and (iq[ncols-1] < iq[ncols-2] or
                          len(m[0]) % ncols != 0):
            ncols = ncols-1

        id = np.reshape(m[0], (-1, ncols)).T[0]
        id = self.__removeTrailingZero(id)
        nrows = len(id)
        if nrows > 1 and id[nrows-1] < id[nrows-2]:
            nrows = nrows-1
        logger.info('psid-psiq ldq %d x %d', nrows, ncols)
        mlen = nrows*ncols
        self.psidq_ldq = {
            'iq': iq[:ncols].tolist(),
            'id': id[:nrows].tolist(),
            'ld': (self.armatureLength*np.reshape(
                m[2][:mlen], (nrows, ncols)).T).tolist(),
            'lq': (self.armatureLength*np.reshape(
                m[3][:mlen], (nrows, ncols)).T).tolist(),
            'psim': (self.armatureLength*np.reshape(
                m[4][:mlen], (nrows, ncols)).T).tolist(),
            'torque': (self.armatureLength*np.reshape(
                m[6][:mlen], (nrows, ncols)).T).tolist()}

    def __read_ldq(self, content):
        "read ld-lq section"
        for i, l in enumerate(content):
            if l.find('[A') > -1:
                break
        m = []
        k = i+2
        for l in content[i+2:]:
            rec = l.split('\t')
            if len(rec) > 7:
                m.append([floatnan(x) for x in rec[:8]])
            elif rec and rec[0].startswith('Curr Id'):
                break
            k += 1

        m = np.array(m).T
        ncols = len(set(m[1]))
        i1 = np.reshape(m[0], (-1, ncols)).T[0]
        nrows = len(i1)
        logger.info('ld-lq %d x %d', nrows, ncols)

        self.ldq = {k: (self.armatureLength*np.reshape(
            v, (nrows, ncols)).T[::-1]).tolist() for k, v in zip(
                ('ld', 'lq', 'psim', 'psid', 'psiq', 'torque'),
                m[2:])}
        self.ldq['beta'] = m[1][:ncols][::-1].tolist()
        self.ldq['i1'] = i1.tolist()

        # skip d-q table
        i = k+3
        for l in content[k+3:]:
            if l.startswith('Losses for'):
                self.__read_losses_tab(content[i:])
            i += 1

    def __read_losses_tab(self, content):
        "read losses of psidq or ldq"
        m = []
        speed = float(content[0].split()[-1])/60.
        logger.info('losses for speed %f', speed)
        nl = 4
        for l in content[4:]:
            rec = l.split('\t')
            if len(rec) == 6:
                m.append([floatnan(x) for x in rec])
            elif rec[0].startswith('P fe'):
                break
            nl += 1
        if not m:
            return
        m = np.array(m).T

        ncols = np.argmax(np.abs(m[1][1:]-m[1][:-1]))+1
        if ncols == 1 and len(m[1]) > 1 and m[1][0] != m[1][1]:
            ncols = 2
        nrows = len(m[2])//ncols
        if ncols * nrows % len(m[3]) != 0:
            if ncols > nrows:
                ncols = ncols-1
            else:
                nrows = nrows-1

        cols = ('styoke', 'stteeth', 'rotor', 'magnet')
        if self.ldq:
            ls = {k: np.reshape(v,
                                (nrows, ncols)).T[::-1].tolist()
                  for k, v in zip(cols, m[2:])}
        else:
            ls = {k: np.reshape(v,
                                (nrows, ncols)).T.tolist()
                  for k, v in zip(cols, m[2:])}
        m = []
        for l in content[nl+3:]:
            rec = l.split('\t')
            if len(rec) == 8:
                m.append([floatnan(x) for x in rec])
            elif not rec and m:
                break

        cols = ('styoke_hyst', 'styoke_eddy',
                'stteeth_hyst', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy')
        if m:
            m = np.array(m).T
            if self.ldq:
                ls.update({k: np.reshape(v,
                                         (nrows, ncols)).T[::-1].tolist()
                           for k, v in zip(cols, m[2:])})
            else:
                ls.update({k: np.reshape(v,
                                         (nrows, ncols)).T.tolist()
                           for k, v in zip(cols, m[2:])})
        ls['speed'] = speed
        if self.ldq:
            self.ldq['losses'] = ls
        elif self.psidq:
            self.psidq['losses'] = ls

    def __read_machine_data(self, content):
        "read machine data section"
        for k in ('beta', 'plfe1', 'plfe2', 'plmag', 'plcu'):
            self.machine[k] = []
        for l in content:
            contentline = l.strip()
            if contentline.startswith('Number of phases'):
                rec = l.split()
                self.machine['m'] = int(rec[-1])
                continue
            if contentline.startswith("Mechanical Power"):
                rec = l.split()
                self.machine['p2'] = floatnan(rec[-1])
                unit = rec[2]
                if unit == '[kW]':
                    self.machine['p2'] *= 1e3
                continue

            for v in [["Phase Current", 'i1'],
                      ["Number of phases", 'm'],
                      ["Current loading", 'A'],
                      ["Current Density", 'J'],
                      ["Cu fillfactor", 'kcu'],
                      ["Therm.Loading", 'AJ'],
                      ["torque", 'torque'],
                      ["Force Density", 'fd'],
                      ['Inductance Ld', 'ld'],
                      ['Inductance Lq', 'lq'],
                      ['Stator Resistance', 'r1'],
                      ['Magn. Flux Psim RMS', 'psim'],
                      ["Speed", 'n'],
                      ['Beta-angle',        'beta'],
                      ["Stator-Fe-Losses noload", 'plfe1_0'],
                      ["Rotor-Fe-Losses noload", 'plfe2_0'],
                      ["Magnet-Losses noload", 'plmag_0'],
                      ["Stator-Fe-Losses  ", 'plfe1'],
                      ["Rotor-Fe-Losses   ", 'plfe2'],
                      ["Magnet-Losses", 'plmag'],
                      ["Cu-Losses", 'plcu'],
                      ["Armature", 'lfe'],
                      ["Efficiency", 'eff']]:

                if contentline.find(v[0]) > -1:
                    si = get_si_factor(contentline)
                    rec = l.split()
                    if v[1] in self.machine and isinstance(
                            self.machine[v[1]], list):
                        self.machine[v[1]].append(si*floatnan(rec[-1]))
                    else:
                        self.machine[v[1]] = si*floatnan(rec[-1])
                    break

        if self.machine['beta'] and len(self.machine['beta']) > 1:
            self.machine['beta'] = self.machine['beta'][1:]
        self.machine['n'] = self.machine['n']/60
        self.machine['lfe'] = 1e-3*self.machine['lfe']
        if self.machine['plfe1']:  # calc sum of losses
            plfe1 = self.machine['plfe1']
            plcu = self.machine.get('plcu', [0.0]*len(plfe1))
            if len(plcu) < len(plfe1):
                self.machine['plcu'] = plcu + [plcu[-1]]*(len(plfe1)-len(plcu))
            self.machine['pltotal'] = [sum(pl)
                                       for pl in zip(*[self.machine[k]
                                                       for k in ('plfe1',
                                                                 'plfe2',
                                                                 'plmag',
                                                                 'plcu')])]
            self.machine['plfe'] = [sum(pl)
                                    for pl in zip(*[self.machine[k]
                                                    for k in ('plfe1',
                                                              'plfe2')])]

    def __read_dq_parameter(self, content):
        if content[1].find('Windings') > -1:  # this is the first section

            for l in content[1:]:
                for v in [['Windings Current', 'i1'],
                          ['Angle I vs. Up', 'beta'],
                          ["LD", 'ld'],
                          ["LQ at nom. current", 'lq'],
                          ["Torque TO", 'torque'],
                          ["Force", 'force'],
                          ["Torque constant Kt", 'kt'],
                          ["Magn.Flux no-load", 'psim0'],
                          ["Voltage Up  (RMS)  no-load", 'up0'],
                          ["Magn.Flux load", 'psim'],
                          ["Voltage Up  load", 'up'],
                          ["Speed", 'speed'],
                          ["Number of Poles", 'npoles'],
                          ["Armature length", 'lfe'],
                          ["Airgap diameter", 'dag']]:
                    if l.find(v[0]) > -1:
                        rec = l.split()
                        self.dqPar[v[1]] = [floatnan(rec[-1])]

            # pdb.set_trace()
            for k in ('speed', 'npoles', 'lfe', 'dag', 'up0', 'up', 'psim0'):
                if k in self.dqPar:
                    self.dqPar[k] = self.dqPar[k][0]
            lfe = self.dqPar['lfe']
            self.dqPar['lfe'] = 1e-3*self.dqPar['lfe']
            for k in ('ld', 'lq', 'psim', 'torque', 'force'):
                if k in self.dqPar:
                    self.dqPar[k][0] = lfe*self.dqPar[k][0]
            if 'torque' in self.dqPar:
                self.dqPar['speed'] = self.dqPar['speed']/60
            if 'dag' in self.dqPar:
                self.dqPar['dag'] = 1e-3*self.dqPar['dag']
            if 'npoles' in self.dqPar:
                self.dqPar['npoles'] = int(self.dqPar['npoles'])
            self.dqPar['i1'] = [self.dqPar['i1'][0]/np.sqrt(2)]
            beta = np.pi*self.dqPar['beta'][0]/180
            iq = np.cos(beta)*self.dqPar['i1'][0]
            id = np.sin(beta)*self.dqPar['i1'][0]
            try:
                w1 = np.pi*self.dqPar['speed']*self.dqPar['npoles']
                r1 = self.machine.get('r1', 0.0)
                uq, ud = (r1*iq + self.dqPar['up'] + id*w1*self.dqPar['ld'][0],
                          r1*id - iq*w1*self.dqPar['lq'][0])
                self.dqPar['u1'] = [np.sqrt(uq**2 + ud**2)]
                self.dqPar['gamma'] = [-np.arctan2(ud, uq)*180/np.pi]
                self.dqPar['psim0'] = lfe*self.dqPar['psim0']
                self.dqPar['phi'] = [self.dqPar['beta'][0] +
                                     self.dqPar['gamma'][0]]
                self.dqPar['cosphi'] = [np.cos(np.pi*phi/180)
                                        for phi in self.dqPar['phi']]
                self.dqPar['i1'].insert(0, 0)
                self.dqPar['u1'].insert(0, self.dqPar.get('up0', 0))
            except KeyError:
                pass

            # if next section is absent
            try:
                self.dqPar['psid'] = [self.dqPar['psim'][0]]
                self.dqPar['psiq'] = [self.dqPar['lq'][0] *
                                      self.dqPar['i1'][-1]]
            except KeyError:
                pass
            return  # end of first section

        # second DQ-Parameter section
        for k in ('i1', 'beta', 'ld', 'lq', 'psim', 'up',
                  'psid', 'psiq', 'torque', 'torque_fe', 'torque_sim',
                  'p2', 'u1_fe', 'u1_sim', 'gamma', 'phi', 'Lho', 'Lh2'):
            self.dqPar[k] = []
        lfe = 1e3*self.dqPar['lfe']

        keys = []
        for l in content:
            rec = self.__findNums(l)
            if len(rec) == 8:
                for k, r in zip(*[keys, rec]):
                    self.dqPar[k].append(floatnan(r))
                self.dqPar['p2'].append(self.dqPar['torque'][-1] *
                                        lfe*2*np.pi*self.dqPar['speed'])
                for k in ('ld', 'lq', 'psim', 'psid', 'up', 'psiq', 'torque'):
                    if self.dqPar[k]:
                        self.dqPar[k][-1] = lfe * self.dqPar[k][-1]
            elif len(rec) == 7:
                self.dqPar['torque_fe'].append(floatnan(rec[2]))
                self.dqPar['torque_sim'].append(floatnan(rec[3]))
                self.dqPar['u1_fe'].append(floatnan(rec[4]))
                self.dqPar['u1_sim'].append(floatnan(rec[5]))
                self.dqPar['gamma'].append(floatnan(rec[6]))
                self.dqPar['phi'].append(floatnan(rec[1]) +
                                         self.dqPar['gamma'][-1])
            elif len(rec) == 5:  # self and mutual inductances
                self.dqPar['Lho'].append(lfe*floatnan(rec[2]))
                self.dqPar['Lh2'].append(lfe*floatnan(rec[3]))
            else:
                headers = l.split()
                try:  # must distinguish ee and pm
                    if headers[5] == 'Voltage':
                        keys = ('i1', 'beta', 'ld', 'lq', 'up',
                                'psid', 'psiq', 'torque')
                    if headers[5] == 'Psi_magn':
                        keys = ('i1', 'beta', 'ld', 'lq', 'psim',
                                'psid', 'psiq', 'torque')
                except:
                    pass
        if self.dqPar['psim']:
            self.dqPar.pop('up', None)
        else:
            self.dqPar.pop('psim', None)
        try:
            w1 = np.pi*self.dqPar['speed']*self.dqPar['npoles']
            r1 = self.machine.get('r1', 0.0)
            beta = np.array(self.dqPar['beta'])/180*np.pi
            iq = np.cos(beta)*self.dqPar['i1']
            id = np.sin(beta)*self.dqPar['i1']
            up = w1*np.array(self.dqPar['psim'])
            ld = np.array(self.dqPar['ld'])
            lq = np.array(self.dqPar['lq'])
            uq = r1*iq + up + id*w1*ld
            ud = r1*id - iq*w1*lq
            self.dqPar['up'] = up.tolist()
            self.dqPar['u1'] = np.sqrt(uq**2 + ud**2).tolist()
            self.dqPar['gamma'] = (-np.arctan2(ud, uq)*180/np.pi).tolist()
            self.dqPar['phi'] = (beta/np.pi*180 + self.dqPar['gamma']).tolist()
            self.dqPar['cosphi'] = [np.cos(np.pi*phi/180)
                                    for phi in self.dqPar['phi']]
            self.dqPar['i1'].insert(0, 0)
            self.dqPar['u1'].insert(0, self.dqPar.get('up0', 0))
        except:
            pass

    def __read_weights(self, content):
        #              Stator-Iron      - Conductors      - Magnets
        #                105.408	     22.542	      0.000
        #              Rotor-Iron       - Conductors       - Magnets
        #                 45.041	      0.000	     17.556
        if self.weights:
            return
        scale = 1  # assume kg unit
        if content[0].split()[-1] == '[gr]':
            scale = 1e-3
        for line in content[2:]:
            rec = line.split()
            if rec[0] != 'Stator-Iron' and rec[0] != 'Rotor-Iron':
                self.weights.append([scale*floatnan(x) for x in rec])

    def __read_areas(self, content):
        #  Area [mm**2]: Stator-Iron           - slots         - Magnets
        #                 16585.2	      2383.1	         0.0
        #  Area [mm**2]: Rotor-Iron            - slots         - Magnets
        #                  6372.7	         0.0	      1930.2
        for line in content:
            rec = line.split()
            if rec and rec[0] != 'Area':
                self.areas.append([floatnan(x) for x in rec])

    def __read_inertia(self, content):
        #  Inertia: GD**2/4 [gr m**2/mm]
        #  Stator :	        0.230195
        #  Rotor  :	        0.011774
        self.inertia = []
        pat = re.compile(r'\[(\w\w) .+\]')
        unit = pat.findall(content[0])
        f = 1
        if unit and unit[-1] == 'gr':
            f = 1e-3
        for line in content[1:]:
            x = line.split()
            if x:
                self.inertia.append(f*floatnan(x[-1]))

    def __read_losses(self, content):
        losses = {}
        # find results for angle:
        for i, l in enumerate(content):
            if l.startswith('Results for Angle I-Up'):
                losses['beta'] = floatnan(l.split(':')[-1])
                losses['current'] = floatnan(
                    content[i+1].split(':')[-1])/np.sqrt(2)
                for k in ('winding', 'staza', 'stajo', 'rotfe',
                          'magnetJ', 'magnetB', 'r1', 'total'):
                    losses[k] = 0.0
                continue

            if l.find('Cu-losses') > -1:
                rec = self.__findNums(content[i+1])
                if len(rec) > 0:
                    losses['winding'] += floatnan(rec[0])
                    losses['total'] += floatnan(rec[0])
                if len(rec) > 2:
                    losses['r1'] += floatnan(rec[2])
                continue

            elif l.startswith('StZa') or l.startswith('RoZa'):
                rec = self.__findNums(content[i+2])
                if len(rec) == 2:
                    losses['stajo'] = floatnan(rec[1])
                    losses['staza'] = floatnan(rec[0])
                    losses['total'] += losses['staza']+losses['stajo']
                continue

            if l.startswith('StJo') or l.startswith('RoJo') or \
               _statloss.search(l):
                rec = self.__findNums(content[i+2])
                if len(rec) == 2:
                    losses['stajo'] = floatnan(rec[0])
                    losses['staza'] = floatnan(rec[1])
                    losses['total'] += losses['staza']+losses['stajo']
                elif len(rec) == 1:
                    t = l.split(':')[-1].strip()
                    if t == 'Iron':
                        losses['staza'] = floatnan(rec[0])
                    else:
                        losses['stajo'] += floatnan(rec[0])
                    losses['total'] += losses['staza']+losses['stajo']
                continue

            if _rotloss.search(l):
                rec = self.__findNums(content[i+2])
                if len(rec) == 1:
                    rotfe = floatnan(rec[0])
                    losses['rotfe'] += rotfe
                    losses['total'] += rotfe
                continue

            if l.find('Fe-Losses-Rotor') > -1:
                rec = self.__findNums(content[i+3])
                if len(rec) == 2:
                    if content[i+1].find('Iron') > -1 and content[i+1].find('StJo') > 0:
                        self.external_rotor = True
                        # TODO: there might be better places to check this
                        losses['stajo'] = floatnan(rec[0])
                        losses['staza'] = floatnan(rec[1])
                        losses['total'] += losses['staza']+losses['stajo']
                    else:
                        losses['rotfe'] = floatnan(rec[1])
                        losses['total'] += losses['rotfe']
                continue

            if l.find('Magnet-Losses') > -1:
                rec = self.__findNums(content[i+1])
                if len(rec) == 1:
                    losses['magnetJ'] = float(rec[0])
                    #losses['magnetB'] = float(Nan)
                if len(rec) == 2:
                    losses['magnetJ'] = float(rec[0])
                    losses['magnetB'] = float(rec[1])
                losses['total'] += losses['magnetJ']

        if 'total' in losses:
            losses['totalfe'] = sum([losses[k] for k in ('staza',
                                                         'stajo',
                                                         'rotfe')])
            self.losses.append(losses)

    def __read_hysteresis_eddy_current_losses(self, content):
        losses = dict(fft=dict(), stator=dict(), rotor=dict())
        part = 'fft'  # either stator or rotor or fft
        k = ''
        for i, l in enumerate(content):
            if l.startswith('*************'):
                for part in losses:
                    self.losses[-1][part] = dict()
                    for k in losses[part]:
                        for x in losses[part][k]:
                            x[0] = int(x[0])
                            if(losses[part][k]):
                                if len(losses[part][k][0]) == 4:
                                    cols = ('order_el', 'freq', 'hyst', 'eddy')
                                else:
                                    cols = ('order_mech', 'order_el',
                                            'freq', 'hyst', 'eddy')

                                self.losses[-1][part][k] = {k1: l
                                                            for k1, l in zip(cols,
                                                                             zip(*losses[part][k]))}
                self.__read_losses(content[i+1:])
                break

            if l.find('StJo') > -1 or \
               l.find('RoJo') > -1:
                k = 'stajo'
            elif l.find('StZa') > -1 or \
                    l.find('StatorIron') > -1 or \
                    l.find('RoZa') > -1:
                k = 'staza'
            elif l.find('Iron') > -1 and l.find('Stator') > -1:
                k = 'staza'
            elif l.find('Iron') > -1:
                if self.external_rotor:
                    k = 'staza'
                else:
                    k = 'rotor'
            elif l.find(': Roto') > -1:
                k = 'rotor'
            elif l.find(': Ring') > -1:
                k = 'ring'
            elif l.find('Stat') > -1 and l.find('Stator: ') < 0:
                k = 'stajo'
            elif l.find('Stator: ') > -1:
                part = 'stator'
                k = l.split(':')[-1].strip()
            elif l.find('Rotor: ') > -1:
                part = 'rotor'
                k = l.split(':')[-1].strip()
            else:
                if k and k not in losses[part]:
                    losses[part][k] = []
                try:
                    rec = self.__findNums(l)
                    if 'nan' in rec:
                        continue
                    if len(rec) == 4:
                        losses[part][k].append([floatnan(x) for x in rec])
                    elif len(rec) == 5:  # FEMAG Rel 8.3 with el/mech order
                        losses[part][k].append([floatnan(x)
                                                for i, x in enumerate(rec)])
                    else:
                        part = 'fft'

                except:
                    pass

    def get(self, name, r=None):
        """return value of key name
        name can be a list such as ['torque[1]', 'ripple']
        or a string: 'dqPar'
        """
        try:
            if isinstance(name, str):  # ignore r
                lname, indx = splitindex(name)
                if indx is None:
                    return self.__getattr__(name)
                return self.__getattr__(lname).__getitem__(indx)

            if len(name) > 1:
                lname, indx = splitindex(name[0])
                if r:
                    if indx is None:
                        return self.get(name[1:],
                                        getattr(r, lname))
                    return self.get(name[1:],
                                    getattr(r, lname).__getitem__(indx))

                if indx is None:
                    return self.get(name[1:],
                                    getattr(self, lname))
                return self.get(name[1:],
                                getattr(self, lname).__getitem__(indx))

            lname, indx = splitindex(name[0])
            if r:
                if indx is None:
                    return r.get(name[0])
                return r.get(lname).__getitem__(indx)

            if indx is None:
                return self.__getattr__(name[0])
            return self.__getattr__(lname).__getitem__(indx)
        except (KeyError, IndexError, AttributeError):
            return None

    def __getattr__(self, k):
        return self.__dict__[k]

    def items(self):
        return [(k, self.get(k)) for k in ('version',
                                           'type',
                                           'filename',
                                           'date',
                                           'areas',
                                           'inertia',
                                           'torque',
                                           'torque_fft',
                                           'psidq',
                                           'psidq_ldq',
                                           'machine',
                                           'lossPar',
                                           'flux',
                                           'flux_fft',
                                           'wdgfactors',
                                           'airgapInduction',
                                           'magnet',
                                           'scData',
                                           'dqPar',
                                           'ldq',
                                           'losses',
                                           'demag',
                                           'linearForce',
                                           'linearForce_fft',
                                           'powerSituation',
                                           'characteristics') if self.get(k)]

    def __str__(self):
        "return string format of this object"
        if self.type:
            return "\n".join([
                'FEMAG {}: {}'.format(self.version, self.type),
                'File: {}  {}'.format(self.filename, self.date)] +
                ['{}: {}'.format(k, v)
                 for k, v in self.items()])
        return "{}"

    def __repr__(self):
        "representation of this object"
        return self.__str__()


def read(filename):
    """Read BCH/BATCH results from file *filename*."""
    import io
    bchresults = Reader()
    with io.open(filename, encoding='latin1', errors='ignore') as f:
        bchresults.read(f.readlines())
    return bchresults


if __name__ == "__main__":
    import json
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    b = read(filename)
    json.dump({k: v for k, v in b.items()}, sys.stdout)
    # print(b)
