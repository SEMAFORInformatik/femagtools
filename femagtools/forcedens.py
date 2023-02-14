# -*- coding: utf-8 -*-
"""
    femagtools.forcedens
    ~~~~~~~~~~~~~~~~~~~~

    Read Force Density Plot Files


"""
import os
import re
import glob
import numpy as np
import logging

logger = logging.getLogger('femagtools.forcedens')

filename_pat = re.compile(r'^(\w+)_(\d{3}).PLT(\d+)')
num_pat = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')
pos_pat = re.compile(r'^\s*POSITION\s*\[(\w+)\]')
unit_pat = re.compile(r'\[([^\]]+)')


def _readSections(f):
    """return list of PLT sections

    sections are surrounded by lines starting with '[***'
    or 2d arrays with 7 columns

    Args:
      param f (file) PLT file to be read

    Returns:
      list of sections
    """

    section = []
    for line in f:
        if line.startswith('[****') or pos_pat.match(line):
            if section:
                if len(section) > 2 and section[1].startswith('Date'):
                    yield section[1:]
                else:
                    yield section
                if line.startswith('[****'):
                    section = []
                else:
                    section = [line.strip()]
        else:
            section.append(line.strip())
    yield section


class ForceDensity(object):

    def __init__(self):
        self.type = ''
        self.positions = []
        pass

    def __read_version(self, content):
        rec = content[0].split(' ')
        if len(rec) > 3:
            self.version = rec[3]
        else:
            self.version = rec[-1]
            
    def __read_project_filename(self, content):
        self.project = content[1].strip()
        
    def __read_nodes_and_mesh(self, content):
        self.nodes, self.elements, self.quality = \
            [float(r[0]) for r in [num_pat.findall(l)
                                   for l in content[:3]]]
        for l in content[3:]:
            m = re.match(r'\*+([^\*]+)\*+', l)
            if m:
                self.type = m.group(1).strip()
                return
            
    def __read_date_and_title(self, content):
        d = content[0].split(':')[1].strip().split()
        dd, MM, yy = d[0].split('.')
        hh, mm = ''.join(d[1:-1]).split('.')
        self.date = '{}-{}-{}T{:02}:{:02}'.format(
            yy, MM, dd, int(hh), int(mm))

        if len(content) > 6:
            self.title = content[2].strip() + ', ' + content[6].strip()
        else:
            self.title = content[2].strip()
            
        self.current = float(num_pat.findall(content[4])[0])
        
    def __read_filename(self, content):
        self.filename = content[0].split(':')[1].strip()
        
    def __read_position(self, content):
        d = dict(position=float(num_pat.findall(content[0])[-1]),
                 unit=unit_pat.findall(content[0].split()[1])[0])
        cols = content[2].split()
        labels = cols[::2] # either X, FN, FT, B_N, B_T
                           # or X FX FY B_X B_Y
        d['column_units'] = {k: u for k, u in zip(labels,
                                                  [unit_pat.findall(u)[0]
                                                   for u in cols[1::2]])}
        m = []
        for l in content[4:]:
            rec = l.split()[1:]
            if len(rec) == len(labels):
                m.append([float(x) for x in rec])
        d.update({k: v for k, v in zip(labels, list(zip(*m)))})

        self.positions.append(d)

    def read(self, filename):
        with open(filename) as f:
            for s in _readSections(f.readlines()):
                logger.debug('Section %s' % s[0:2])
                if s[0].startswith('FEMAG'):
                    self.__read_version(s)
                elif s[0].startswith('Project'):
                    self.__read_project_filename(s)
                elif s[0].startswith('Number'):
                    self.__read_nodes_and_mesh(s)
                elif s[0].startswith('File'):
                    self.__read_filename(s)
                elif s[0].startswith('Date'):
                    self.__read_date_and_title(s)
                elif s[0].startswith('POSITION'):
                    self.__read_position(s)

    def _prepare(self): 
        """return FN and FT Matrix"""
        ntiles = int(360/self.positions[0]['X'][-1])
        try:
            FN = np.tile(
                np.array([p['FN'][:-1] for p in self.positions[:-1]]).T,
                (ntiles, 1))
            FT = np.tile(
                np.array([p['FT'][:-1] for p in self.positions[:-1]]).T,
                (ntiles, 1))
        except AttributeError:
            return []
        return [FN, FT]
    
    def fft(self):
        """return 2D FFT of the Force Density"""
        fdens = self._prepare()
        FN, FT= fdens[0], fdens[1] 
        N = FN.size
        dim = FN.shape[0]//2
        # Dimension [Position °] X [Time Steps s]
        FN_MAG = np.fft.fft2(FN)[0:dim, :]/N
        FT_MAG = np.fft.fft2(FT)[0:dim, :]/N

        return dict(fn_harm=dict(amplitude=np.abs(FN_MAG), 
                                phase=np.angle(FN_MAG)), 
                    ft_harm=dict(amplitude=np.abs(FT_MAG), 
                                phase=np.angle(FT_MAG)) 
                    )
    
    def items(self):
        return [(k, getattr(self, k)) for k in ('version',
                                                'type',
                                                'title',
                                                'current',
                                                'filename',
                                                'date',
                                                'positions')]

    def __str__(self):
        "return string format of this object"
        if self.type:
            return "\n".join([
                'FEMAG {}: {}'.format(self.version, self.type),
                'File: {}  {}'.format(self.filename, self.date)] +
                             ['{}: {}'.format(k, v)
                              for k, v in self.items()])
        return "{}"


def read(filename):
    f = ForceDensity()
    f.read(filename)
    return f


def readall(workdir='.'):
    """collect all recent PLT files
    returns list of ForceDensity objects
    """
    plt = dict()
    pltfiles = sorted(glob.glob(os.path.join(workdir, '*_*.PLT*')))
    base = os.path.basename(pltfiles[-1])
    lastserie = filename_pat.match(base).groups()[1]
    for p in pltfiles:
        base = os.path.basename(p)
        m = filename_pat.match(base)
        if m and lastserie == m.groups()[1]:
            model, i, k = m.groups()
            fdens = ForceDensity()
            fdens.read(p)
            logging.info("%s: %s", p, fdens.title)
            if model in plt:
                plt[model].append(fdens)
            else:
                plt[model] = [fdens]
    return plt


if __name__ == "__main__":
    import matplotlib.pyplot as pl
    import femagtools.plot
    import sys
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    fdens = read(filename)

    # Show the results

    title = '{}, Rotor position {}'.format(
        fdens.title, fdens.positions[0]['position'])
    pos = fdens.positions[0]['X']
    FT_FN = (fdens.positions[0]['FT'],
             fdens.positions[0]['FN'])
    femagtools.plot.forcedens(title, pos, FT_FN)    
    pl.show()

    title = 'Force Density Harmonics'
    femagtools.plot.forcedens_fft(title, fdens)
    pl.show()
     
