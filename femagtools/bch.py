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
        

class Reader:
    """Reads a BCH/BATCH-File"""
    _numPattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)')

    def __init__(self):
        self._fft = None
        self.type = None
        self.filename = ''
        self.project = ''
        self.date = ''
        self.version = ''
        # known types:
        # MULTIPLE CALCULATION OF FORCES AND FLUX
        # Fast cogging calculation OF FORCES AND FLUX
        # Fast LD-LQ-Identification
        # Fast Psid-Psiq-Identification
        # Fast PM-Synchronous-Motor Simulation
        self.wdg = None
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
        self.scData = {}
        self.dqPar = {}
        self.ldq = {}
        self.losses = []
        self.demag = {}
        self.weights = []
        self.weight = {}
        self.areas = []
        self.current_angles = []
        self.dispatch = {
            'General Machine Data': Reader.__read_general_machine_data,
            'Weigths': Reader.__read_weights,
            'Number of Nodes': Reader.__read_nodes_and_mesh,
            'Windings Data': Reader.__read_dummy,
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
            'Losses [W]':  Reader.__read_losses,
            'Losses for speed [1/min]': Reader.__read_losses_tab,
            'Losses from PSID-Psiq-Identification for speed [1/min]':
            Reader.__read_losses_tab,
            'Magnet loss data': Reader.__read_dummy,
            'Project File name': Reader.__read_project_filename,
            'File name': Reader.__read_filename,
            'Windings input data': Reader.__read_windings,
            'Control parameters for Loss calculation': Reader.__read_lossPar,
            'PSID-Psiq-Identification': Reader.__read_psidq,
            'Ld-Lq-Identifikation aus PSID-Psiq-Identification': Reader.__read_psidq_ldq,
            'Ld-Lq-Identification RMS-values': Reader.__read_ldq,
            'Machine Data Rotor': Reader.__read_dummy,
            'Current Angles defined from no-load test': Reader.__read_current_angles,
            'FEMAG Version': Reader.__read_version,
            'Simulation Data': Reader.__read_dummy,
            'Area [mm**2]': Reader.__read_areas,
            'Basic Machine parameters': Reader.__read_dummy,
            'Winding': Reader.__read_dummy,
            'Calculation time [sec]': Reader.__read_dummy,
            'Results for Angle I-Up [Degree]': Reader.__read_dummy,
            'Demagnetisation': Reader.__read_demagnetization,
            'Transient short circuit': Reader.__read_short_circuit,
            'Flux observed': Reader.__read_flux}

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
            title = s[0].split(':')[0].strip()
            if title == 'Function':
                title = s[0].split(':')[1].strip()
            logger.debug("'%s': %d", title, len(s[1:]))
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
                    wdgs.append((int(rec[1]), floatnan(rec[-1])))
            else:
                index = False
        self.windings = []
        for m in set([int(abs(w)) for w in np.array(wdgs).T[0]]):
            self.windings.append([w for w in wdgs if abs(w[0]) == m])

    def __read_dummy(self, content):
        return
    
    def __read_current_angles(self, content):
        self.current_angles = []
        for l in content:
            rec = self._numPattern.findall(l)
            if len(rec) == 3:
                self.current_angles.append(floatnan(rec[-1]))
        return
    
    def __read_lossPar( self,content ):
        self.lossPar={
            'fo':[],
            'Bo':[],
            'ch':[],
            'cw':[],
            'hf':[],
            'ef':[],
            'ic':[],
            'gamfe':[],
            'thetaw':[],
            'fillfe':[],
            'mat':[]}
        for l in content:
            for v in [ ["Base Frequency", 'fo' ],
                       ["Base Induction", 'Bo' ],
                       ["Hysteresis-Coefficient", 'ch'],
                       ["Eddycurrent-Coefficient", 'cw'],
                       ["Hysteresis-Frequency-Coefficient", 'hf'],
                       ["Eddycurrent-Frequency-Coefficient", 'ef'],
                       ["Induction-Coefficient", 'ic'],
                       ["Specific Weight Iron", 'gamfe'],
                       ["Conductor Temperature", 'thetaw'],
                       ["Fillfactor Iron", 'fillfe'],
                       ["Material factor", 'mat'] ]:

                if l.find(v[0])>-1:
                    rec=self._numPattern.findall(l)
                    if len(rec)>0:
                        self.lossPar[v[1]].append(floatnan(rec[-1]))
                    break

        return l
    
    def __read_demagnetization( self, content ):
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

    def __read_short_circuit(self, content):
        "read short circuit section"
        if content[2].startswith('Time'):
            m=[]
            for l in content:
                rec = l.split()
                if len(rec) == 5 and not rec[0].startswith('Time'):
                    m.append([floatnan(x) for x in rec])

            m = np.array(m).T
            self.scData['time'] = m[0]
            self.scData['ia'] = m[1]
            self.scData['ib'] = m[2]
            self.scData['ic'] = m[3]
            self.scData['torque'] = m[4]
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
                
    def __read_flux(self, content):
        "read and append flux section"

        f = {'displ': [], 'flux_k': [], 'voltage_dpsi': [],
             'voltage_four': [], 'current_k': [], 'voltage_ir': []}

        for l in content:
            rec = l.split()
            if len(rec) == 7:
                f['displ'].append(floatnan(rec[1].strip()))
                f['flux_k'].append(floatnan(rec[2].strip()))
                f['voltage_dpsi'].append(floatnan(rec[3].strip()))
                f['voltage_four'].append(floatnan(rec[4].strip()))
                f['current_k'].append(floatnan(rec[5].strip()))
                f['voltage_ir'].append(floatnan(rec[6].strip()))
            elif len(rec) > 0 and rec[0].startswith('['):
                f['displunit'] = re.search(r"\[([^\]]*)\]", l).group(1).strip()
            elif l.startswith('Flux-Area'):
                self.wdg = rec[-1]
                if self.wdg not in self.flux:
                    self.flux[self.wdg] = []
        self.flux[self.wdg].append(f)
        self._fft = Reader.__read_flux_fft
        
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
        if ncols > 1 and iq[ncols-1] < iq[ncols-2]:
            ncols = ncols-1

        id = np.reshape(m[0], (-1, ncols)).T[0]
        nrows = len(id)
        if nrows > 1 and id[nrows-1] < id[nrows-2]:
            nrows = nrows-1
        logger.info('psid-psiq %d x %d', nrows, ncols)
        mlen = nrows*ncols
        self.psidq = {
            'iq': iq[:ncols].tolist(),
            'id': id[:nrows].tolist(),
            'psid': (self.armatureLength*np.reshape(
                m[3][:mlen],
                (nrows, ncols))).T.tolist(),
            'psiq': (self.armatureLength*np.reshape(
                m[4][:mlen],
                (nrows, ncols))).T.tolist(),
            'torque': (self.armatureLength*np.reshape(
                m[6][:mlen],
                (nrows, ncols))).T.tolist()}
        
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
        if ncols > 1 and iq[ncols-1] < iq[ncols-2]:
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
        for l in content[i+2:]:
            rec = l.split('\t')
            if len(rec) > 7:
                m.append([floatnan(x) for x in rec[:8]])
                
        m = np.array(m).T
        ncols = len(set(m[1]))
        i1 = np.reshape(m[0], (-1, ncols)).T[0]
        nrows = len(i1)
        logger.info('ld-lq %d x %d', nrows, ncols)

        self.ldq = {
            'beta': m[1][:ncols][::-1].tolist(),
            'i1': i1.tolist(),
            'ld': (self.armatureLength*np.reshape(
                m[2], (nrows, ncols)).T[::-1]).tolist(),
            'lq': (self.armatureLength*np.reshape(
                m[3], (nrows, ncols)).T[::-1]).tolist(),
            'psim': (self.armatureLength*np.reshape(
                m[4], (nrows, ncols)).T[::-1]).tolist(),
            'psid': (self.armatureLength*np.reshape(
                m[5], (nrows, ncols)).T[::-1]).tolist(),
            'psiq': (self.armatureLength*np.reshape(
                m[6], (nrows, ncols)).T[::-1]).tolist(),
            'torque': (self.armatureLength*np.reshape(
                m[7], (nrows, ncols)).T[::-1]).tolist()}

    def __read_losses_tab(self, content):
        "read losses of psidq or ldq"
        m = []
        for l in content[4:]:
            rec = l.split('\t')
            if len(rec) == 6:
                m.append([floatnan(x) for x in rec])
        m = np.array(m).T
        ncols = len(set(m[1]))
        nrows = len(m[2])//ncols
        l = dict(
            styoke=np.reshape(m[2],
                              (nrows, ncols)).T[::-1].tolist(),
            stteeth=np.reshape(m[3],
                               (nrows, ncols)).T[::-1].tolist(),
            rotor=np.reshape(m[4],
                             (nrows, ncols)).T[::-1].tolist(),
            magnet=np.reshape(m[5],
                              (nrows, ncols)).T[::-1].tolist())
        if self.ldq:
            self.ldq['losses'] = l
        else:
            self.psidq['losses'] = l
 
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
                    if v[1] in self.machine and isinstance(
                            self.machine[v[1]], list):
                        self.machine[v[1]].append(si*floatnan(rec[-1]))
                    else:
                        self.machine[v[1]] = si*floatnan(rec[-1])
                    break

        if 'beta' in self.machine:
            self.machine['beta'] = self.machine['beta'][1:]
        self.machine['n'] = self.machine['n']/60
        self.machine['lfe'] = 1e-3*self.machine['lfe']
        i1 = self.machine['i1']
        if 'plfe1' in self.machine:  # calc sum of losses
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
        if content[1].find('Windings') > -1:
            for l in content[1:]:
                for v in [['Windings Current', 'i1'],
                          ['Angle I vs. Up', 'beta'],
                          ["LD", 'ld'],
                          ["LQ at nom. current", 'lq'],
                          ["Torque TO", 'torque'],
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

            for k in ('speed', 'npoles', 'lfe', 'dag', 'up0', 'up', 'psim0'):
                if k in self.dqPar:
                    self.dqPar[k] = self.dqPar[k][0]
            lfe = self.dqPar['lfe']
            self.dqPar['lfe'] = 1e-3*self.dqPar['lfe']
            self.dqPar['dag'] = 1e-3*self.dqPar['dag']
            for k in ('ld', 'lq', 'psim', 'torque'):
                self.dqPar[k][0] = lfe*self.dqPar[k][0]
            self.dqPar['speed'] = self.dqPar['speed']/60
            self.dqPar['npoles'] = int(self.dqPar['npoles'])
            self.dqPar['i1'] = [self.dqPar['i1'][0]/np.sqrt(2)]
            beta = np.pi*self.dqPar['beta'][0]/180
            iq = np.cos(beta)*self.dqPar['i1'][0]
            id = np.sin(beta)*self.dqPar['i1'][0]
            w1 = np.pi*self.dqPar['speed']*self.dqPar['npoles']
            uq, ud = (self.dqPar['up'] + id*w1*self.dqPar['ld'][0],
                      iq*w1*self.dqPar['lq'][0])
            self.dqPar['u1'] = [np.sqrt(uq**2 + ud**2)]
            self.dqPar['gamma'] = [-np.arctan2(ud, uq)*180/np.pi]
            self.dqPar['psim0'] = lfe*self.dqPar['psim0']
            self.dqPar['phi'] = [self.dqPar['beta'][0] +
                                 self.dqPar['gamma'][0]]
            self.dqPar['cosphi'] = [np.cos(np.pi*phi/180)
                                    for phi in self.dqPar['phi']]
            self.dqPar['i1'].insert(0, 0)
            self.dqPar['u1'].insert(0, self.dqPar['up0'])
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
                losses['current'] = floatnan(content[i+1].split(':')[-1])
                losses['total'] = 0.0
                for k in ('stacu', 'staza', 'stajo', 'rotfe',
                          'magnetJ', 'magnetB'):
                    losses[k] = 0.0
                break

        for i, l in enumerate(content):
            if l.find('Cu-losses') > -1:
                rec = self._numPattern.findall(content[i+1])
                if len(rec) > 0:
                    losses['winding'] = floatnan(rec[0])
                    losses['total'] += losses['winding']
                if len(rec) > 2:
                    losses['r1'] = floatnan(rec[2])
                    
            if l.startswith('StZa') or l.startswith('RoZa'):
                rec = self._numPattern.findall(content[i+2])
                if len(rec) == 2:
                    losses['stajo'] = floatnan(rec[1])
                    losses['staza'] = floatnan(rec[0])
                    losses['total'] += losses['staza']+losses['stajo']

            if l.startswith('StJo') or l.startswith('RoJo'):
                rec = self._numPattern.findall(content[i+2])
                if len(rec) == 2:
                    losses['stajo'] = floatnan(rec[0])
                    losses['staza'] = floatnan(rec[1])
                    losses['total'] += losses['staza']+losses['stajo']
                        
            if l.find('Fe-Losses   Rotor:') > -1 or \
               l.find('Fe-Losses   Stator:') > -1:
                rec = self._numPattern.findall(content[i+2])
                if len(rec) == 1:
                    losses['rotfe'] = floatnan(rec[0])
                    losses['total'] += losses['rotfe']
                    
            if l.find('Fe-Losses-Rotor') > -1:
                rec = self._numPattern.findall(content[i+3])
                if len(rec) == 2:
                    losses['rotfe'] = floatnan(rec[1])
                    losses['total'] += losses['rotfe']
                    
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
                self.__read_losses(content[i:])
                break

            if l.find('StJo') > -1 or \
               l.find('RoJo') > -1 or \
               l.find('Stat') > -1:
                k = 'stajo'
            elif l.find('StZa') > -1 or \
                 l.find('RoZa') > -1:
                k = 'staza'
            elif l.find('Iron') > -1 or \
                 l.find('Roto') > -1:
                k = 'rotor'
            else:
                try:
                    rec = self._numPattern.findall(l)
                    if len(rec) == 4:
                        losses[k].append([floatnan(x) for x in rec])
                except:
                    pass
                
        for k in losses:
            for x in losses[k]:
                x[0] = int(x[0])
        self.losses[-1]['fft'] = losses
    
    def get(self, name, r=None):
        "return value of key name"

        if isinstance(name, str):
            return self.__getattr__(name)
        if r and type(r) == dict:
            for k in name:
                r = r.get(k)
            return r
        if len(name) > 1:
            if r:
                self.get(name[1:], getattr(r, name[0]))
            return self.get(name[1:], getattr(self, name[0]))
        return self.__getattr__(name[0])
        
    def __getattr__(self, k):
        return self.__dict__[k]

    def items(self):
        return [
            ('version', self.version),
            ('type', self.type),
            ('filename', self.filename),
            ('date', self.date),
            ('torque', self.torque),
            ('torque_fft', self.torque_fft),
            ('psidq', self.psidq),
            ('psidq_ldq', self.psidq_ldq),
            ('machine', self.machine),
            ('lossPar', self.lossPar),
            ('flux', self.flux),
            ('flux_fft', self.flux_fft),
            ('airgapInduction', self.airgapInduction),
            ('magnet', self.magnet),
            ('scData', self.scData),
            ('dqPar', self.dqPar),
            ('ldq', self.ldq),
            ('losses', self.losses),
            ('demag', self.demag)]

    def __str__(self):
        "return string format of this object"
        if self.type:
            return "\n".join([
                'FEMAG {}: {}'.format(self.version, self.type),
                'File: {}  {}'.format(self.filename, self.date),
                'torque:{}'.format(self.torque),
                'torque_fft:{}'.format(self.torque_fft),
                'psidq: {}'.format(self.psidq),
                'psidq_ldq: {}'.format(self.psidq_ldq),
                'machine: {}'.format(self.machine),
                'lossPar: {}'.format(self.lossPar),
                'flux: {}'.format(self.flux),
                'flux_fft: {}'.format(self.flux_fft),
                'magnet: {}'.format(self.magnet),
                'airgapInduction: {}'.format(self.airgapInduction),
                'scData: {}'.format(self.scData),
                'dqPar: {}'.format(self.dqPar),
                'ldq: {}'.format(self.ldq),
                'losses: {}'.format(self.losses),
                'demag: {}'.format(self.demag)])

        return "{}"
    
    def __repr__(self):
        "representation of this object"
        return self.__str__()

def main():
    from io import open
#    with open('logging.json', 'rt') as f:
#        logging.config.dictConfig( json.load(f) )
    bch = Reader()
    for name in sys.argv[1:]:
        with codecs.open(name, encoding='ascii') as f:
            bch.read( f )
        print( bch.type )
        print( bch.date )
        #print bch.torque
        print( bch.losses )
        #print( bch.losses[-1]['stajo'] + bch.losses[-1]['stajo'] )
        #print( bch.areas )
        #print( bch.weights )
        #print( bch.windings )
        #print bch.psidq['id']
        #print bch.psidq['iq']
        #print( bch.machine )
        #print( bch.dqPar['beta'] )
        #print( bch.dqPar['ld'] )
        #print( bch.dqPar )
        #print( bch.psidq )
        #print( bch.ldq )
        #print( bch.flux['1'][0] )#[0]['current_k'] )
        #print( bch.flux_fft['1'][0] )#[0]['current_k'] )
        #print( bch.torque_fft )
        #print( bch.scData['time'] )
        #d={}
        #bch.get(['weight','magnet'])
        #for k in bch.psidq:
        #    d[k]=bch.psidq[k].tolist()
            
        #json.dump(d, sys.stdout)
        #print bch.getStep()
        plot=False
        if plot:
            import matplotlib.pyplot as pl
            for k in ('1','2','3'):
                pl.plot( bch.flux[k][0]['displ'], bch.flux[k][0]['current_k'] )
            pl.xlabel('Displ. / Deg')
            pl.ylabel('Current / A')
            pl.grid()
            pl.show()
        
    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
