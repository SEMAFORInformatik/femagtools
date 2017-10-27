# -*- coding: utf-8 -*-
"""
    femagtools.poc
    ~~~~~~~~~~~~~~

    Manage POC files



"""


class Poc:
    def __init__(self, arg, parameters=dict()):
        """initialize this object from a set of parameters or a file

        Args:
            arg filename or pole_pitch
        """
        for k in parameters.keys():
            self.__setattr__(k, parameters[k])

        if isinstance(arg, str):
            with open(arg) as f:
                self.readfile(f)
        else:
            self.pole_pitch = arg
            self.pocType = parameters.get('pocType', 'Function')
            self.shape_current = parameters.get('shape_current', 'sin')
            self.num_winding = parameters.get('num_winding', 3)
            self.key_winding = parameters.get('key_winding',
                                          list(range(1, self.num_winding+1)))
            b = parameters.get('offset', 0)
            self.phi_voltage_winding = parameters.get('phi_voltage_winding',
                                                      [b+i*360/self.num_winding
                                                       for i in range(self.num_winding)])

    def __setattr__(self, name, val):
        self.__dict__[name] = val  # this will create the attribute name

    def write(self, pocfilename):
        """create a new pocfile and write the data"""
        pocfile = open(pocfilename, mode='w')

        self.writefile(pocfile)

        pocfile.close()

    def writefile(self, pocfile):
        # windings (num and keys)
        pocfile.write("{0}\n".format(self.num_winding))

        for idx in self.key_winding[:self.num_winding]:
            pocfile.write("{0}\n".format(idx))

        # phi voltage windings
        for v in self.phi_voltage_winding[:self.num_winding]:    
            pocfile.write("{0}\n".format(v))

        # rest
        if self.pole_pitch:
            pocfile.write("{0}\n".format(self.pole_pitch))
        if self.pocType in ['fun', 'har', 'hsp']:
            pocfile.write("{0}\n{1}\n".format(
                self.pocType, self.func_steps))
            for i, val in enumerate(self.func_current[:self.func_steps]):
                if self.harmonic_id[i] and \
                   self.pocType == 'hsp':
                    pocfile.write("{0}, {1}, {2}".format(
                        self.harmonic_id[i],
                        val, self.func_phi[i]))
                pocfile.write("{0}, {1}\n".format(
                    self.func_phi[i], val))

        if self.pocType == 'Function':
            pocfile.write("{0}\n".format(self.shape_current))

        if 'skew_angle' in self.__dict__:
            pocfile.write("{0}\n".format(self.skew_angle))
        if 'num_skew_steps' in self.__dict__:
            pocfile.write("{0}\n".format(self.num_skew_steps))
        pocfile.write("\n")

    def readfile(self, pocfile):
        """read poc file"""
        self.num_winding = int(pocfile.readline())
        self.key_winding = []
        for i in range(self.num_winding):
            self.key_winding.append(int(pocfile.readline()))
        self.phi_voltage_winding = []
        for v in self.key_winding:
            self.phi_voltage_winding.append(float(pocfile.readline()))
        self.pole_pitch = float(pocfile.readline())
        self.pocType = pocfile.readline().strip()
        if self.pocType in ['fun', 'har', 'hsp']:
            self.func_current = []
            self.func_phi = []
            self.func_steps = int(pocfile.readline())
            if self.pocType == 'hsp':
                self.harmonic_id = []
            for i in range(self.func_steps):
                l = pocfile.readline().strip().split(',')
                if len(l) > 2:
                    self.harmonic_id.append(int(l[0]))
                    self.func_current.append(float(l[1]))
                    self.func_phi.append(float(l[2]))
                else:
                    self.func_current.append(float(l[0]))
                    self.func_phi.append(float(l[1]))
        else:
            self.shape_current=self.pocType
            self.pocType='Function'
        
        try:
            self.skew_angle=float(pocfile.readline())
            self.num_skew_steps=int(pocfile.readline())
        except ValueError:
            pass

    def getProps( self ):
        keys=['num_winding',
              'key_winding',
              'phi_voltage_winding',
              'pole_pitch',
              'pocType',
              'shape_current',
              'skew_angle',
              'num_skew_steps',
              'func_steps',
              'harmonic_id',
              'func_phi',
              'func_current']
        props={}
        for k in keys:
            if k in self.__dict__:
                props[k]=self.__dict__[k]
        return props

def curr(x, n, A, phi ):
    "return fourier sum"
    if isinstance(A,list):
        amax=max(A)
        s=np.zeros(len(x))
        for ai,ni, phii in zip(A, n, phi):
            #if abs(ai/amax)>1e-2:
            s += curr(x, ni, ai, phii)
        return s
    return A*np.sin(n*x-phi)

if __name__ == "__main__":
    p = Poc('2p_sin.poc')
    print(p.getProps())

    p = Poc(60)
    print(p.getProps())

#    import numpy as np
#    import numpy.linalg as la
#    import matplotlib.pyplot as pl

#    A = [0.9866, 0.011181, 0.018624, 0.022322, 0.020109, 0.016582]
#    n = [1,11,13,14,16,17]
#    phi = [0, 225.801, 162.3245, 195.1087, 321.5113, -3.517]
#    x = np.linspace( 0, 2*np.pi, 40 )
#    y = curr( x, n, A,phi)
#    pl.plot( x, y )
#    pl.grid()
#    pl.show()
