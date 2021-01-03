# -*- coding: utf-8 -*-
"""
    femagtools.windings
    ~~~~~~~~~~~~~~~~~~~

    Handling windings

 Conventions

 Number of slots: Q
 Numper of pole pairs: p
 Number of phases: m
 Number of layers: l
 Number of wires per slot side: n
 Number of parallel circuits: g
 Number of slots per pole and phase: q = Q/p/2/m
 Number of windings per phase: w1 = Q * n * l/2/m/g
"""
import numpy as np
import femagtools.bch


class Windings(object):

    def __init__(self, arg):
        if isinstance(arg, femagtools.bch.Reader):
            self.m = arg.machine['m']
            self.Q = arg.machine['Q']
            self.p = arg.machine['p']
            self.windings = arg.windings
        else:
            for k in arg.keys():
                setattr(self, k, arg[k])

    def sequence(self):
        """returns sequence of winding keys"""
        return list(zip(*sorted([(k, self.windings[k]['PHI'][0])
                                 for k in self.windings.keys()],
                                key=lambda wdg: wdg[1])))[0]

    def slots(self, key):
        """returns slot indexes of winding key"""
        ngen = self.m*self.Q//np.gcd(self.Q, self.m*2*self.p)
        taus = 360/self.Q
        startpos = int(self.windings[key]['PHI'][0]/taus)
        s = [int(x/taus) for x in self.windings[key]['PHI']]
        layers = 1 if len(s) == len(set(s)) else 2

        dim = int(layers*ngen/self.m)
        slots = [int(x/taus) + ngen*n for n in range(self.Q//ngen)
                 for x in self.windings[key]['PHI'][:dim]]
        return np.array(slots).reshape((np.gcd(self.Q, self.p), -1)) - startpos

    def axis(self):
        """returns axis angle of winding 1 in mechanical system"""
        return self.current_linkage()['alfa0']
    
    def current_linkage(self, k=1):
        taus = 2*np.pi/self.Q
        t = np.gcd(self.Q, self.p)
        slots = self.slots(k)[0]
        dirs = self.windings[k]['dir']
        curr = np.concatenate([np.array(dirs)*(1 - 2*(n % 2)) 
                               for n in range(len(slots)//len(dirs))])

        NY=4096
        y = np.zeros(NY*self.Q//t)
        for i in range(self.Q//t):
            if i in set(slots):
                y[NY*i+NY//2] = np.sum(curr[slots==i])
        yy = [np.sum(y[:i+1]) for i in range(0, len(y))]                
        yy[:NY//2] = yy[-NY//2:]
        yy = np.tile(yy-np.mean(yy),t)
        yy /= np.max(yy)
        #y = np.tile(y,t)

        N = len(yy)
        Y = np.fft.fft(yy)
        i = np.argmax(np.abs(Y[:N//2]))
        a = 2*np.abs(Y[i])/N
        freq = np.fft.fftfreq(N, d=taus/NY)
        T0 = np.abs(1/freq[i])
        alfa0 = np.angle(Y[i])
        #if alfa0 < 0: alfa0 += 2*np.pi
        pos_fft = np.linspace(0, self.Q/t*taus)
        D = (a*np.cos(2*np.pi*pos_fft/T0+alfa0))
        return dict(
            pos = [i*taus/NY for i in range(len(y))],
            current_linkage=yy[:NY*self.Q//t].tolist(),
            alfa0=-alfa0/self.p,
            pos_fft = pos_fft.tolist(),
            current_linkage_fft=D.tolist())

        
if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if sys.argv[1:]:
        bch = femagtools.bch.read(sys.argv[1])
        wdgs = Windings(bch)
    else:
        testdata=[
            dict(Q = 90, p = 12, m = 3,
                 windings = {1: {
                     'dir': [-1, 1, 1, -1, -1, -1, 1, 1, -1, -1],
                     'N': [8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
                     'PHI': [2.0, 6.0, 6.0, 18.0, 18.0, 22.0, 34.0, 34.0, 50.0, 50.0]}}),

            dict(Q = 54, p = 6, m = 3,
                 windings = {1: {
                     'dir': [1, 1, 1, -1, -1, -1],
                     'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                     'PHI': [3.3333, 3.3333, 10.0, 30.0, 36.6666, 36.6666]}}),

            dict(Q = 168, p = 7, m = 3,
                 windings = {1: {'dir': [1, 1, 1, 1, 1, 1, -1, -1],
                                 'N': [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
                                 'PHI': [1.0714, 1.0714, 3.2143, 3.2143, 5.3572, 7.5, 22.5001, 24.6429]},
                             2: {'dir': [1, 1, 1, 1, 1, 1, 1, 1],
                                 'N': [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
                                 'PHI': [13.9286, 16.0715, 18.2143,  18.2143,  20.3572, 20.3572, 22.5001, 24.6429]},
                             3: {'dir': [-1, -1, -1, -1, -1, -1, -1, -1],
                                 'N': [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
                                 'PHI': [5.3572, 7.5, 9.6429, 9.6429, 11.7857, 11.7857, 13.9286, 16.0715]}})]
        wdgs = Windings(testdata[0])

    c = wdgs.current_linkage()
    #print('alfa0={0:6.3f}'.format(wdgs.axis()/np.pi*180))

    plt.title('Q={0}, p={1}, alfa0={2:6.3f}'.format(wdgs.Q, wdgs.p, c['alfa0']/np.pi*180))
    plt.plot(np.array(c['pos'])/np.pi*180, c['current_linkage'])
    plt.plot(np.array(c['pos_fft'])/np.pi*180, c['current_linkage_fft'])

    phi = [c['alfa0']/np.pi*180, c['alfa0']/np.pi*180]
    y = [min(c['current_linkage_fft']), 1.1*max(c['current_linkage_fft'])]
    plt.plot(phi, y, '--')
    plt.annotate("", xy=(phi[0], y[0]), 
                 xytext=(0, y[0]), arrowprops=dict(arrowstyle="->"))

    plt.grid()
    plt.show()
