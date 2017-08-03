# -*- coding: utf-8 -*-
"""
    femagtools.bch
    ~~~~~~~~~~~~~~

    Reading BCH/BATCH files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np
import re
import codecs
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
    pattern = re.compile("\[([A-Za-z/0-9]+)\]")
    search = pattern.search(contentline)
    if search:
        if search.group(1).startswith('kW'):
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
                yield section
                section = []
        else:
            section.append(line.strip())
    yield section


class Reader:
    """Reads a BCH/BATCH-File"""
    _numPattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')

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
        self.scData = {}
        self.dqPar = {}
        self.ldq = {}
        self.losses = []
        self.demag = {}
        self.weights = []
        self.weight = {}
        self.characteristics = {}
        self.areas = []
        self.current_angles = []
        self.demagnetization = {}
        self.dispatch = {
            'General Machine Data': Reader.__read_general_machine_data,
            'Weigths': Reader.__read_weights,
            'Number of Nodes': Reader.__read_nodes_and_mesh,
            'Windings Data': Reader.__read_dummy,
            'Winding-Factors': Reader.__read_winding_factors,
            'Machine Data': Reader.__read_machine_data,
            'CAD-Parameter Data': Reader.__read_dummy,
            'Torque-Force': Reader.__read_torque_force,
            'Fourier Analysis': Reader.__read_fft,
            'Machine excitation': Reader.__read_dummy,
            'DQ-Parameter for open Winding Modell': Reader.__read_dq_parameter,
            'Magnet Data': Reader.__read_magnet_data,
            'Date': Reader.__read_date,
            'Inertia': Reader.__read_dummy,
            'Losses': Reader.__read_dummy,
            'Fe-Hysteresis- and Eddy current Losses[W] > 0.1 % Max':
            Reader.__read_hysteresis_eddy_current_losses,
            'Losses [W]': Reader.__read_losses,
            'Losses for speed [1/min]': Reader.__read_losses_tab,
            'Losses from PSID-Psiq-Identification for speed [1/min]':
            Reader.__read_losses_tab,
            'Magnet loss data': Reader.__read_dummy,
            'Project File name': Reader.__read_project_filename,
            'File name': Reader.__read_filename,
            'Windings input data': Reader.__read_windings,
            'Control parameters for Loss calculation': Reader.__read_lossPar,
            'PSID-Psiq-Identification': Reader.__read_psidq,
            'Ld-Lq-Identifikation aus PSID-Psiq-Identification':
            Reader.__read_psidq_ldq,
            'Ld-Lq-Identification RMS-values': Reader.__read_ldq,
            'Machine Data Rotor': Reader.__read_dummy,
            'Current Angles defined from no-load test':
            Reader.__read_current_angles,
            'FEMAG Version': Reader.__read_version,
            'Simulation Data': Reader.__read_simulation_data,
            'Area [mm**2]': Reader.__read_areas,
            'Basic Machine parameters': Reader.__read_dummy,
            'Winding': Reader.__read_dummy,
            'Calculation time [sec]': Reader.__read_dummy,
            'Results for Angle I-Up [Degree]': Reader.__read_dummy,
            'Demagnetisation': Reader.__read_demagnetization,
            'Transient short circuit': Reader.__read_short_circuit,
            'Flux observed': Reader.__read_flux,
            'Linear Force': Reader.__read_linear_force,
            'Demagnetization Data': Reader.__read_demagnetization,
            '** Characteristics of Permanent-Magnet-Motors **':
            Reader.__read_characteristics}

    def getStep(self):
        """@returns displacement step of flux values"""
        if len(self.flux) > 0:
            return self.flux[0]['displ'][1]-self.flux[0]['displ'][0]
        return None

    def read(self, lines):
        """read bch file

        Args:
          lines (list of str) the text lines of the BCH file
        """
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

            if title in self.dispatch:
                self.dispatch[title](self, s)

        if len(self.weights) > 0:
            w = list(zip(*self.weights))
            self.weight['iron'] = sum(w[0])
            self.weight['conductor'] = sum(w[1])
            self.weight['magnet'] = sum(w[2])
            self.weight['total'] = sum([sum(l) for l in w])

    def __read_version(self, content):
        self.version = content[0].split(' ')[3]

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
            [floatnan(r[0]) for r in [self._numPattern.findall(l)
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
                    rec = self._numPattern.findall(l)
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
                rec = self._numPattern.findall(line)
                if len(rec) == 4:
                    self.wdgfactors[-1]['order'].append(int(rec[0]))
                    self.wdgfactors[-1]['wfac'].append(float(rec[1]))
                    self.wdgfactors[-1]['skewf'].append(float(rec[2]))
                    self.wdgfactors[-1]['total'].append(float(rec[3]))
                    
    def __read_dummy(self, content):
        return

    def __read_simulation_data(self, content):
        for line in content:
            if line.startswith('Number of Phases m'):
                self.machine['m'] = int(float((line.split()[-1])))
        return

    def __read_current_angles(self, content):
        self.current_angles = []
        for l in content:
            rec = self._numPattern.findall(l)
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
                      ["Material factor", 'mat']]:

                if l.find(v[0]) > -1:
                    rec = self._numPattern.findall(l)
                    if len(rec) > 0:
                        self.lossPar[v[1]].append(floatnan(rec[-1]))
                    break

        return l

    def __read_demagnetization(self, content):
        keys = ('displ', 'current_1', 'current_2', 'current_3',
                'h_max', 'h_av', 'area')
        
        for l in content:
            rec = self._numPattern.findall(l)
            if len(rec) == 7:
                for r, k in zip(rec, keys):
                    if k in self.demag:
                        self.demag[k].append(floatnan(r.strip()))
                    else:
                        self.demag[k] = [floatnan(r.strip())]
            else:
                for v in [["Limit Hc value", "lim_hc"],
                          ["Max. Magnetization", "br_max"],
                          ["Min. Magnetization", "br_min"],
                          ["Area demagnetized", "area"]]:
                    if l.find(v[0]) > -1:
                        rec = self._numPattern.findall(l)
                        if len(rec) > 0:
                            self.demag[v[1]] = floatnan(rec[-1])

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
                
    def __read_characteristics(self, content):
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
                    rec = self._numPattern.findall(l)
                    if len(rec) > 0:
                        self.characteristics[v[1]] = floatnan(rec[-1])
                        break

        self.characteristics['ldq'] = {}
        self.armatureLength = self.characteristics['lfe']
        m = []
        for k, l in enumerate(content[i+3:]):
            if l.startswith('[[***'):
                break
            rec = l.split('\t')
            if len(rec) == 6:
                m.append([floatnan(x) for x in rec])

        if m:
            m = np.array(m).T
            ncols = len(set(m[1]))
            i1 = np.reshape(m[0], (-1, ncols)).T[0]
            nrows = len(i1)

            logger.info('characteristics ld-lq %d x %d', nrows, ncols)
            self.characteristics['ldq'] = {
                'beta': m[1][:ncols][::-1].tolist(),
                'i1': i1.tolist(),
                'ld': (self.armatureLength*np.reshape(
                    m[2], (nrows, ncols)).T[::-1]).tolist(),
                'lq': (self.armatureLength*np.reshape(
                    m[3], (nrows, ncols)).T[::-1]).tolist(),
                'psim': (self.armatureLength*np.reshape(
                    m[4], (nrows, ncols)).T[::-1]).tolist(),
                'torque': (self.armatureLength*np.reshape(
                    m[5], (nrows, ncols)).T[::-1]).tolist()}

        m = []
        columns = [['n', 'id', 'iq', 'torque', 'p2'],
                   ['beta', 'cos_phi', 'u1', 'um'],
                   ['lang', 'ud', 'uq', 'i1'],
                   ['lang', 'ld', 'lq', 'psim']]
        nsec = 0
        self.characteristics['speed_torque'] = {}
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
                    self.characteristics['speed_torque'][k] = m[j].tolist()
                m = []
                nsec += 1
            else:
                rec = self._numPattern.findall(l)
                if len(rec) > 3:
                    m.append([floatnan(x) for x in rec])

    def __read_flux(self, content):
        "read and append flux section"

        f = {'displ': [], 'flux_k': [], 'voltage_dpsi': [],
             'voltage_four': [], 'current_k': [], 'voltage_ir': []}

        for l in content:
            rec = l.split()
            if l.startswith('Flux-Area'):
                areas = self._numPattern.findall(l)
                if not areas:
                    continue

                self.wdg = areas[0] if len(areas)==1 else '{}-{}'.format(
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
        f = {'displ': [], 'magnet_1': [], 'force_x': [],
             'force_y': [], 'f_idpsi': []}

        for l in content:
            rec = self._numPattern.findall(l)
            if len(rec) > 4:
                f['displ'].append(floatnan(rec[1].strip()))
                f['magnet_1'].append(floatnan(rec[2].strip()))
                f['force_x'].append(floatnan(rec[3].strip()))
                f['force_y'].append(floatnan(rec[4].strip()))
                # TODO f['f_idpsi'].append(floatnan(rec[5].strip()))
        if len(f['displ']) > 0:
            ripple_x = max(f['force_x']) - min(f['force_x'])
            ripple_y = max(f['force_y']) - min(f['force_y'])
            f['ripple_x'] = ripple_x
            f['ripple_y'] = ripple_y
            self.linearForce.append(f)
        
        self._fft = Reader.__read_linearForce_fft

    def __read_linearForce_fft(self, content):
        "read and append linear force fft section"
        if not self._fft:
            return
        linearForce_fft = dict(order=[], force=[], force_perc=[],
                               a=[], b=[])
        for l in content:
            rec = self._numPattern.findall(l)
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
            rec = self._numPattern.findall(l)
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
            rec = self._numPattern.findall(l)
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
                rec = self._numPattern.findall(l)
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
        ncols = len(set(m[1]))
        iq = m[1]
        if ncols > 1 and (iq[ncols-1] < iq[ncols-2] or
                          len(m[0]) % ncols != 0):
            ncols = ncols-1

        id = np.reshape(m[0], (-1, ncols)).T[0]
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
        ncols = len(set(m[1]))
        iq = m[1]
        if ncols > 1 and (iq[ncols-1] < iq[ncols-2] or
                          len(m[0]) % ncols != 0):
            ncols = ncols-1

        id = np.reshape(m[0], (-1, ncols)).T[0]
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
        for l in content[4:]:
            rec = l.split('\t')
            if len(rec) == 6:
                m.append([floatnan(x) for x in rec])
        if not m:
            return
        m = np.array(m).T
        try:
            ncols = len(set(m[1]))
            nrows = len(m[2])//ncols
            if ncols * nrows % len(m[3]) != 0:
                if ncols > nrows:
                    ncols = ncols-1
                else:
                    nrows = nrows-1

            l = {k: np.reshape(v,
                               (nrows, ncols)).T[::-1].tolist()
                 for k, v in zip(('styoke', 'stteeth', 'rotor', 'magnet'),
                                 m[2:])}
            l['speed'] = speed
            if self.ldq:
                self.ldq['losses'] = l
            else:
                self.psidq['losses'] = l
        except:
            pass
        
    def __read_machine_data(self, content):
        "read machine data section"
        for k in ('beta', 'plfe1', 'plfe2', 'plmag'):
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
                    # Check if torque is saved in kNm instead of Nm. thomas.maier/OSWALD
                    if v[1] == "torque" and rec[1] == "[kNm]":
                        si = si * 1000
                    if v[1] in self.machine and isinstance(
                            self.machine[v[1]], list):
                        self.machine[v[1]].append(si*floatnan(rec[-1]))
                    else:
                        self.machine[v[1]] = si*floatnan(rec[-1])
                    break

        if self.machine['beta']:
            self.machine['beta'] = self.machine['beta'][1:]
        self.machine['n'] = self.machine['n']/60
        self.machine['lfe'] = 1e-3*self.machine['lfe']
        i1 = self.machine['i1']
        if self.machine['plfe1']:  # calc sum of losses
            self.machine['i1'] = i1*len(self.machine['plfe1'])
            plfe1 = self.machine['plfe1']
            plcu = self.machine.get('plcu', 0.0)
            self.machine['plcu'] = [plcu]*len(plfe1)
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
        import pdb
        if content[1].find('Windings') > -1:
            
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
                self.dqPar['u1'].insert(0, self.dqPar['up0'])
            except KeyError:
                pass

            # if next section is absent
            try:
                self.dqPar['psid'] = [self.dqPar['psim'][0] * np.sqrt(2.)]
                self.dqPar['psiq'] = [self.dqPar['lq'][0] * self.dqPar['i1'][-1]
                                      * np.sqrt(2.)]
            except KeyError:
                pass
            return
        
        for k in ('i1', 'beta', 'ld', 'lq', 'psim', 'psid', 'psiq', 'torque',
                  'p2', 'u1', 'gamma', 'phi'):
            self.dqPar[k] = []
        lfe = 1e3*self.dqPar['lfe']

        for l in content:
            rec = self._numPattern.findall(l)
            if len(rec) == 8:
                for k, r in zip(*[('i1', 'beta', 'ld', 'lq', 'psim',
                                   'psid', 'psiq', 'torque'), rec]):
                    self.dqPar[k].append(floatnan(r))
                self.dqPar['p2'].append(self.dqPar['torque'][-1] *
                                        lfe*2*np.pi*self.dqPar['speed'])
                for k in ('ld', 'lq', 'psim', 'psid', 'psiq', 'torque'):
                    self.dqPar[k][-1] = lfe * self.dqPar[k][-1]
            elif len(rec) == 7:
                self.dqPar['u1'].append(floatnan(rec[4]))
                self.dqPar['gamma'].append(floatnan(rec[6]))
                self.dqPar['phi'].append(self.dqPar['beta'][-1] +
                                         self.dqPar['gamma'][-1])

        self.dqPar['cosphi'] = [np.cos(np.pi*phi/180)
                                for phi in self.dqPar['phi']]
        self.dqPar['i1'].insert(0, 0)
        self.dqPar['u1'].insert(0, self.dqPar['up0'])
        
    def __read_weights( self, content ):
    #              Stator-Iron      - Conductors      - Magnets 
    #                105.408	     22.542	      0.000
    #              Rotor-Iron       - Conductors       - Magnets  
    #                 45.041	      0.000	     17.556
        if self.weights:
            return
        for line in content[2:]:
            rec=line.split()
            if rec[0] != 'Stator-Iron' and rec[0] != 'Rotor-Iron':
                self.weights.append([floatnan(x) for x in rec] )

    def __read_areas( self, content ):
    #Area [mm**2]: Stator-Iron           - slots         - Magnets  
    #                 16585.2	      2383.1	         0.0
    #Area [mm**2]: Rotor-Iron            - slots         - Magnets 
    #                  6372.7	         0.0	      1930.2
        for line in content:
            rec=line.split()
            if rec[0] != 'Area':
                self.areas.append([floatnan(x) for x in rec] )

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
                rec = self._numPattern.findall(content[i+1])
                if len(rec) > 0:
                    losses['winding'] += floatnan(rec[0])
                    losses['total'] += floatnan(rec[0])
                if len(rec) > 2:
                    losses['r1'] += floatnan(rec[2])
                continue
                    
            elif l.startswith('StZa') or l.startswith('RoZa'):
                rec = self._numPattern.findall(content[i+2])
                if len(rec) == 2:
                    losses['stajo'] = floatnan(rec[1])
                    losses['staza'] = floatnan(rec[0])
                    losses['total'] += losses['staza']+losses['stajo']
                continue
                
            if l.startswith('StJo') or l.startswith('RoJo') or \
               _statloss.search(l):
                rec = self._numPattern.findall(content[i+2])
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
                rec = self._numPattern.findall(content[i+2])
                if len(rec) == 1:
                    losses['rotfe'] = floatnan(rec[0])
                    losses['total'] += losses['rotfe']
                continue
                    
            if l.find('Fe-Losses-Rotor') > -1:
                rec = self._numPattern.findall(content[i+3])
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
                rec = self._numPattern.findall(content[i+1])
                if len(rec) == 1:
                    losses['magnetJ'] = float(rec[0])
                    #losses['magnetB'] = float(Nan)
                if len(rec) == 2:
                    losses['magnetJ'] = float(rec[0])
                    losses['magnetB'] = float(rec[1])
                losses['total'] += losses['magnetJ']
                
        if 'total' in losses:
            self.losses.append(losses)

    def __read_hysteresis_eddy_current_losses(self, content):
        losses = dict(staza=[], stajo=[], rotor=[])
        for i, l in enumerate(content):
            if l.startswith('*************'):
                for k in losses:
                    for x in losses[k]:
                        x[0] = int(x[0])
                self.losses[-1]['fft'] = {k: {k1: l
                                              for k1, l in zip(['order',
                                                                'freq',
                                                                'hyst',
                                                                'eddy'],
                                                               zip(*losses[k]))}
                                          for k in losses}
                self.__read_losses(content[i+1:])
                break

            if l.find('StJo') > -1 or \
               l.find('RoJo') > -1:
                k = 'stajo'
            elif l.find('StZa') > -1 or \
                 l.find('RoZa') > -1:
                k = 'staza'
            elif l.find('Iron') > -1 or \
                 l.find('Roto') > -1:
                if self.external_rotor:
                    k = 'staza'
                else:
                    k = 'rotor'
            elif l.find('Stat') > -1:
                k = 'stajo'
            else:
                try:
                    rec = self._numPattern.findall(l)
                    if len(rec) == 4:
                        losses[k].append([floatnan(x) for x in rec])
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

def main():
#    from io import open
#    with open('logging.json', 'rt') as f:
#        logging.config.dictConfig( json.load(f) )
    bch = Reader()
    for name in sys.argv[1:]:
        with codecs.open(name, encoding='ascii') as f:
            bch.read(f)
        print(bch.type)
        print(bch.date)
        print(bch.characteristics)
        # print(bch.torque)
        # print(bch.losses[-1])
        # print(bch.linearForce)
        # print( bch.losses[-1]['stajo'] + bch.losses[-1]['stajo'] )
        # print( bch.areas )
        # print( bch.weights )
        # print( bch.windings )
        # print( bch.psidq['id'] )
        # print( bch.psidq['iq'] )
        # print( bch.psidq_ldq['psim'] )
        # print( bch.machine )
        # print( bch.dqPar['beta'] )
        # print( bch.dqPar['ld'] )
        # print( bch.dqPar )
        # print( bch.psidq )
        # print( bch.ldq )
        # print( bch.flux['1'][0] )#[0]['current_k'] )
        # print( bch.flux_fft['1'][0] )#[0]['current_k'] )
        # print( bch.torque_fft )
        # print( bch.scData['time'] )
        # d={}
        # bch.get(['weight','magnet'])
        # for k in bch.psidq:
        #    d[k]=bch.psidq[k].tolist()

        # json.dump(d, sys.stdout)
        # print bch.getStep()
        plot=False
        if plot:
            import matplotlib.pyplot as pl
            for k in ('1', '2', '3'):
                pl.plot( bch.flux[k][0]['displ'], bch.flux[k][0]['current_k'] )
            pl.xlabel('Displ. / Deg')
            pl.ylabel('Current / A')
            pl.grid()
            pl.show()

    return 0

if __name__ == "__main__":
    import json
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    b = Reader()

    with codecs.open(filename, encoding='ascii') as f:
        b.read(f)

    #json.dump(b, sys.stdout)
    print(b)
#    status = main()
#    sys.exit(status)
