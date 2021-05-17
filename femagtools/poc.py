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
        self.skew_angle = 0.0
        self.num_skew_steps = 0
        for k in parameters.keys():
            self.__setattr__(k, parameters[k])

        if isinstance(arg, str):
            with open(arg) as f:
                self.readfile(f)
        else:
            self.pole_pitch = arg
            self.pocType = parameters.get('pocType', 'Function')
            self.shape_current = parameters.get('shape_current', 'sin')
            self.key_winding = parameters.get('key_winding',
                                              [str(i+1) for i in range(3)])
            b = parameters.get('offset', 0)
            num_winding = len(self.key_winding)
            self.phi_voltage_winding = parameters.get('phi_voltage_winding',
                                                      [b+i*360/num_winding
                                                       for i in range(num_winding)])

    def __setattr__(self, name, val):
        self.__dict__[name] = val  # this will create the attribute name

    def filename(self):
        prefix = self.pocType
        if prefix == 'Function':
            prefix = self.shape_current
        return '{}_{}p.poc'.format(
            prefix, round(720/self.pole_pitch))

    def write(self, pocfilename):
        """create a new pocfile and write the data"""
        with open(pocfilename, mode='w') as pocfile:
            pocfile.write('\n'.join(self.content()))

    def content(self):
        """return list of lines"""
        # windings (num and keys)
        num_winding = len(self.key_winding)
        content = ["{0}".format(num_winding)]
        content.extend(["{0}".format(k)
                        for k in self.key_winding[:num_winding]])

        # phi voltage windings
        content.extend(["{0}".format(v)
                       for v in self.phi_voltage_winding[:num_winding]])

        # rest
        if self.pole_pitch:
            content.append("{0}".format(self.pole_pitch))
        if self.pocType in ['fun', 'har', 'hsp']:
            func_steps = len(self.func_current)
            content += [f"{self.pocType}", f"{func_steps}"]
            if (self.pocType == 'fun' and
                    num_winding*func_steps > len(self.func_current)):
                self.func_current = num_winding * self.func_current
                self.func_phi = num_winding * self.func_phi
                func_steps = num_winding*func_steps
            for i, val in enumerate(self.func_current[:func_steps]):
                if self.pocType == 'hsp':
                    content.append("{0}, {1}, {2}".format(
                        self.harmonic_id[i],
                        val, self.func_phi[i]))
                elif self.pocType == 'har':
                    content.append("{0}, {1}".format(
                        val, self.func_phi[i]))
                elif self.pocType == 'fun':
                    content.append("{0}, {1}".format(
                        self.func_phi[i], val))

        if self.pocType == 'Function':
            content.append("{0}".format(self.shape_current))

        content += [f"{self.skew_angle}",
                    f"{int(self.num_skew_steps)}",
                    '']
        return content

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
            func_steps = int(pocfile.readline())
            if self.pocType == 'hsp':
                self.harmonic_id = []
            for i in range(func_steps):
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
        keys=['key_winding',
              'phi_voltage_winding',
              'pole_pitch',
              'pocType',
              'shape_current',
              'skew_angle',
              'num_skew_steps',
              'harmonic_id',
              'func_phi',
              'func_current']
        props={}
        for k in keys:
            if k in self.__dict__:
                props[k]=self.__dict__[k]
        return props

class HspPoc(Poc):
    def __init__(self, harm=[1], amp=[1], phi=[0]):
        Poc.__init__(self, 0, dict(
            pocType='hsp', harmonic_id=harm,
            func_current=amp, func_phi=phi))

        
def curr_har(x, n, A, phi ):
    "return fourier sum"
    import numpy as np
    if isinstance(A, list):
        return np.sum([curr_har(x, *nAp)
                       for nAp in zip(n, A, phi)], axis=0)
    return A*np.sin(n*x - phi)


def curr_fun(phi, A, phi_wdg ):
    "return curr"
    import numpy as np
    if np.asarray(A).shape == 1:
        return [[np.asarray(phi)+x, A] for x in phi_wdg] 
    return [phi, A]


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
