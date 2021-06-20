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
 Number of slots per pole and phase: q = Q/p/2/m
 Number of coils per phase: c = Q * l/2/m
 Number of parallel circuits (coil groups): g
 Number of windings per phase: w1 = Q * n * l/2/m/g
"""
import numpy as np
import femagtools.bch


def q1q2yk(Q, p, m, l=1):
    """returns q1, q2, Yk, k1, k2"""
    if l == 1:  # single layer
        t = np.gcd(Q//2, p)
    else:
        t = np.gcd(Q, p)
    Qb = Q//t
    pb = p//t
    if Qb % 2:  # odd
        q2 = (Qb + m)//(2*m) - 1
        q1 = q2 + 1
    else:
        q2 = (Qb)//(2*m)
        q1 = q2
    if l == 2:
        n = 1
        while (n*Qb + 1) % pb:
            n += 1
        Yk = (n*Qb + 1)//pb
    else:
        n = 1
        while (n*Qb//2 + 1) % pb:
            n += 1
        Yk = (n*Qb//2 + 1)//pb

    k1 = [(q1 + q2)*i for i in range(m)]
    k2 = (q1*(m+1) + q2*(m-1))//2
    return q1, q2, Yk, k1, k2, Qb

#         pos = [sorted([Yk*(k + n) % Qb
#                       for n in range(q1)]) for k in k1]
#         neg = [sorted([Yk*(k + n + k2) % Qb
#                       for n in range(q2)]) for k in k1]


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

        if hasattr(self, 'windings'):
            return

        l = 1
        if hasattr(self, 'l'):
            l = self.l
        coilwidth = 0
        if l > 1:
            coilwidth = self.Q//self.p//2

        if hasattr(self, 'coilwidth'):
            coilwidth = self.coilwidth

        q1, q2, Yk, k1, k2, Qb = q1q2yk(self.Q, self.p, self.m, l)
        pos = [sorted([(Yk*(k + n)) % Qb for n in range(q1)]) for k in k1]
        neg = [sorted([Yk*(k + n + k2) % Qb for n in range(q2)]) for k in k1]
        taus = 360/self.Q
        taup = self.Q//self.p//2
        if coilwidth:
            slots = [[(k, 1, 1) for k in s] for s in pos]
            for i, s in enumerate(neg):
                slots[i] = sorted(slots[i] + [(k, -1, 1)
                                              for k in s], key=lambda x: x[0])
        else:
            slots = [[(k % (Qb//2), -1 if k % (Qb//2) < i*Qb//3 else 1, 1)
                      for k in s] for i, s in enumerate(pos)]
            for i, s in enumerate(slots):
                slots[i] = sorted(s + [((k[0]+taup) % Qb, -k[1], 1)
                                       for k in s], key=lambda x: x[0])
        if coilwidth:
            slots = [sorted(s + [((k[0]+coilwidth) % (Qb), -1*k[1], 0)
                                 for k in s], key=lambda x:x[0]) for s in slots]
        self.windings = {i+1:  dict(dir=[k[1] for k in s],
                                    N=[1]*len(s), R=[k[2] for k in s],
                                    PHI=[taus/2+k[0]*taus for k in s])
                         for i, s in enumerate(slots)}

    def sequence(self):
        """returns sequence of winding keys"""
        return list(zip(*sorted([(k, self.windings[k]['PHI'][0])
                                 for k in self.windings.keys()],
                                key=lambda wdg: wdg[1])))[0]

    def slots(self, key):
        """returns slot indexes of winding key"""
        ngen = self.m*self.Q//np.gcd(self.Q, self.m*2*self.p)
        taus = 360/self.Q
        s = [round((x-taus/2)/taus)
             for x in self.windings[key]['PHI']]
        layers = 1 if len(s) == len(set(s)) else 2

        dim = int(layers*ngen/self.m)
        slots = [round((x-taus/2)/taus) + 1 + ngen*n
                 for n in range(self.Q//ngen)
                 for x in self.windings[key]['PHI'][:dim]]
        return np.array(slots).reshape((np.gcd(self.Q, self.p), -1))

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

        NY = 4096
        y = np.zeros(NY*self.Q//t)
        for i in range(1, self.Q//t+1):
            if i in set(slots):
                y[NY*(i-1)+NY//2] = np.sum(curr[slots == i])
        yy = [np.sum(y[:i+1]) for i in range(0, len(y))]
        yy[:NY//2] = yy[-NY//2:]
        yy = np.tile(yy-np.mean(yy), t)
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
        pos_fft = np.linspace(0, self.Q/t*taus, self.p//t*60)
        D = (a*np.cos(2*np.pi*pos_fft/T0+alfa0))
        return dict(
            pos=[i*taus/NY for i in range(len(y))],
            current_linkage=yy[:NY*self.Q//t].tolist(),
            alfa0=-alfa0/self.p,
            pos_fft=pos_fft.tolist(),
            current_linkage_fft=D.tolist())

    def zoneplan(self):
        taus = 360/self.Q
        dphi = 1e-3
        q = self.Q/self.p/self.m/2
        t = np.gcd(self.Q, self.p)
        Qb = self.Q//t
        slots = {k: [round((x-taus/2)/taus)
                     for x in self.windings[k]['PHI']]
                 for k in self.windings}
        layers = 1
        avgr = 0
        maxr, minr = max(self.windings[1]['R']), min(self.windings[1]['R'])
        if maxr-minr > 1e-6:
            layers = 2
            avgr = (maxr+minr)/2

            def is_upper(r, phi):
                return r > avgr
        elif len(slots[1]) > len(set(slots[1])):
            layers = 2

            def is_upper(r, phi):
                return phi < -dphi
        else:
            def is_upper(r, phi):
                return True

        upper = [[s+1 for s, x, r in zip(
            slots[key],
            self.windings[key]['PHI'],
            self.windings[key]['R'])
            if is_upper(r, s*taus - (x-taus/2))]
            for key in self.windings]
        udirs = [[d for s, d, x, r in zip(
            slots[key],
            self.windings[key]['dir'],
            self.windings[key]['PHI'],
            self.windings[key]['R'])
            if is_upper(r, s*taus - (x-taus/2))]
            for key in self.windings]
        lower = []
        ldirs = []
        if layers > 1:
            lower = [[s+1 for s, x, r in zip(
                slots[key],
                self.windings[key]['PHI'],
                self.windings[key]['R'])
                if not is_upper(r, s*taus - (x-taus/2))]
                for key in self.windings]
            ldirs = [[d for s, d, x, r in zip(
                slots[key],
                self.windings[key]['dir'],
                self.windings[key]['PHI'],
                self.windings[key]['R'])
                if not is_upper(r, s*taus - (x-taus/2))]
                for key in self.windings]
        if len([n for l in upper for n in l]) < Qb:
            upper = [m + [n+Qb//2 for n in m] for m in upper]
            lower = [m + [n+Qb//2 for n in m] for m in lower]
            udirs = [m + [-d for d in m] for m in udirs]
            ldirs = [m + [-d for d in m] for m in ldirs]

        return ([[d*s for s, d in zip(u, ud)] for u, ud in zip(upper, udirs)],
                [[d*s for s, d in zip(l, ld)] for l, ld in zip(lower, ldirs)])


if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if sys.argv[1:]:
        bch = femagtools.bch.read(sys.argv[1])
        wdgs = Windings(bch)
    else:
        testdata = [
            dict(Q=90, p=12, m=3,
                 windings={1: {
                     'dir': [-1, 1, 1, -1, -1, -1, 1, 1, -1, -1],
                     'N': [8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
                     'PHI': [2.0, 6.0, 6.0, 18.0, 18.0, 22.0, 34.0, 34.0, 50.0, 50.0]}}),

            dict(Q=54, p=6, m=3,
                 windings={1: {
                     'dir': [1, 1, 1, -1, -1, -1],
                     'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                     'PHI': [3.3333, 3.3333, 10.0, 30.0, 36.6666, 36.6666]}}),

            dict(Q=168, p=7, m=3,
                 windings={1: {'dir': [1, 1, 1, 1, 1, 1, -1, -1],
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
    # print('alfa0={0:6.3f}'.format(wdgs.axis()/np.pi*180))

    plt.title('Q={0}, p={1}, alfa0={2:6.3f}'.format(
        wdgs.Q, wdgs.p, c['alfa0']/np.pi*180))
    plt.plot(np.array(c['pos'])/np.pi*180, c['current_linkage'])
    plt.plot(np.array(c['pos_fft'])/np.pi*180, c['current_linkage_fft'])

    phi = [c['alfa0']/np.pi*180, c['alfa0']/np.pi*180]
    y = [min(c['current_linkage_fft']), 1.1*max(c['current_linkage_fft'])]
    plt.plot(phi, y, '--')
    plt.annotate("", xy=(phi[0], y[0]),
                 xytext=(0, y[0]), arrowprops=dict(arrowstyle="->"))

    plt.grid()
    plt.show()
