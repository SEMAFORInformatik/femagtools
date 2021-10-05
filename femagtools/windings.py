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
 Number of turns per phase: w1 = Q * n * l/2/m/g
"""
import numpy as np
import femagtools.bch
from xml.etree import ElementTree as ET

coil_color = ['lime', 'gold', 'magenta',
              'blue', 'brown', 'blueviolet']


def num_basic_windings(Q, p, l):
    """return number of basic windings"""
    if l == 1:  # single layer
        return np.gcd(Q//2, p)
    return np.gcd(Q, p)


def q1q2yk(Q, p, m, l=1):
    """returns q1, q2, Yk, Qb"""
    t = num_basic_windings(Q, p, l)
    Qb = Q//t
    qqb = Qb if l == 2 else Qb//2
    pb = p//t
    if qqb//m % 2:  # odd
        q2 = (qqb + m)//(2*m) - 1
        q1 = q2 + 1
    else:
        q2 = (qqb)//(2*m)
        q1 = q2
    n = 1
    while (n*qqb + 1) % pb:
        n += 1
    Yk = (n*qqb + 1)//pb
    return q1, q2, Yk, Qb


class Winding(object):

    def __init__(self, arg):
        """create winding either from bch winding section or winding dict()

        """
        if isinstance(arg, femagtools.bch.Reader):
            self.m = arg.machine['m']
            self.Q = arg.machine['Q']
            self.p = arg.machine['p']
            self.windings = arg.windings
        else:
            for k in arg.keys():
                setattr(self, k, arg[k])

        self.q = self.Q/2/self.p/self.m  # number of coils per pole and phase

        if hasattr(self, 'windings'):
            # calculate coil width yd and num layers l
            taus = 360/self.Q
            if 'slots' in self.windings[1]:  # custom winding def
                for k in self.windings:
                    w = self.windings[k]
                    w['dir'] = [1 if s > 0 else -1 for s in w['slots']]
                    w['PHI'] = [(2*abs(s)-1)*taus/2 for s in w['slots']]
                    w['R'] = [0 if l == 1 else 1 for l in w['layer']]
            k = 0
            for w in self.windings:
                try:
                    k = self.windings[w]['dir'].index(
                        -self.windings[w]['dir'][0])
                    wk = w
                except ValueError:
                    pass
            if k == 0:
                self.yd = max(self.Q//self.p//2, 1)
            else:
                self.yd = round((self.windings[wk]['PHI'][k] -
                                 self.windings[wk]['PHI'][0])/taus)

            slots = [round((x-taus/2)/taus)
                     for x in self.windings[1]['PHI']]
            self.l = 2
            if len(slots) == len(set(slots)):
                self.l = 1
            return

        layers = 1
        if hasattr(self, 'l'):
            layers = self.l
        else:
            self.l = layers
        coilwidth = max(self.Q//self.p//2, 1)
        if hasattr(self, 'coilwidth'):   # obsolete, use yd instead
            coilwidth = self.coilwidth
        elif hasattr(self, 'yd'):
            coilwidth = self.yd

        self.yd = coilwidth

        q1, q2, Yk, Qb = q1q2yk(self.Q, self.p, self.m, self.l)
        k1 = [(q1 + q2)*i for i in range(self.m)]
        k2 = (q1*(self.m+1) + q2*(self.m-1))//2
        j = 2 if layers == 1 else 1
        pos = [[(j*Yk*(k + n)) % Qb
                for n in range(q1)] for k in k1]
        neg = [[j*Yk*(k + n + k2) % Qb
                for n in range(q2)] for k in k1]
        if layers > 1:
            slots = [sorted([(k, 1, 1)
                             for k in p] + [(k, -1, 1)
                                            for k in n])
                     for n, p in zip(neg, pos)]
            for i, p in enumerate(slots):
                slots[i] = sorted(slots[i] +
                                  [((k[0]+coilwidth) % Qb, -k[1], 0)
                                   for k in slots[i]], key=lambda s: s[0])
        else:
            if (coilwidth + 1) % 2:
                coilwidth += 1
            xneg = [sorted([s for s in n if s+1 % 2] +
                           [(s + coilwidth) % Qb for s in p if s+1 % 2])
                    for n, p in zip(neg, pos)]
            xpos = [sorted([s for s in p if s+1 % 2] +
                           [(s + coilwidth) % Qb for s in n if s+1 % 2])
                    for n, p in zip(neg, pos)]

            slots = [sorted([(k, 1, 1) for k in p] + [(k, -1, 1) for k in n])
                     for n, p in zip(xneg, xpos)]

        if self.m > 3:  # TODO: check this hack
            slots = slots[:self.m//2] + [[((k[0]+1) % self.Q, k[1], k[2])
                                          for k in s]
                                         for s in slots[:self.m//2]]
        taus = 360/self.Q
        self.windings = {
            i+1:  dict(dir=[k[1] for k in s],
                       N=[1]*len(s), R=[k[2] for k in s],
                       PHI=[taus/2+k[0]*taus for k in s])
            for i, s in enumerate(slots)}

    def kw_order(self, n):
        """return winding factor harmonics"""
        if n == 0:
            return self.p
        g = np.arange(-n, n, 1)
        t = num_basic_windings(self.Q, self.p, self.l)
        return self.p + g * self.m*t

    def kwp(self, n=0):
        """pitch factor"""
        nue = n if n and not np.isscalar(n) else self.kw_order(n)
        return np.sin(nue*self.yd*np.pi/self.Q)

    def kwd(self, n=0):
        """zone (distribution) factor"""
        q1, q2, Yk, Qb = q1q2yk(self.Q, self.p, self.m, self.l)
        nue = n if n and not np.isscalar(n) else self.kw_order(n)
        if q1 == q2:
            x = nue*np.pi/self.Q
            return np.sin(q1*x)/(q1*np.sin(x))
        x = nue*np.pi*Yk/self.Q
        k = 2 if self.l == 1 else 1
        return -((np.sin(k*x*q1) - np.cos(x*Qb)*np.sin(k*x*q2)) /
                 ((q1+q2)*np.sin(k*x)))

    def kw(self, n=0):
        """return winding factor"""
        return self.kwp(n) * self.kwd(n)

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

        dim = int(self.l*ngen/self.m)
        slots = [(round((x-taus/2)/taus) + ngen*n) % self.Q + 1
                 for n in range(self.Q//ngen)
                 for x in self.windings[key]['PHI'][:dim]]
        return np.array(slots).reshape((np.gcd(self.Q, self.p), -1))

    def axis(self):
        """returns axis angle of winding 1 in mechanical system"""
        return self.mmf()['alfa0']

    def mmf(self, k=1):
        """returns the dimensionless magnetomotive force (ampere-turns/turns/ampere) and
        winding angle of phase k (rad)"""
        taus = 2*np.pi/self.Q
        t = np.gcd(self.Q, self.p)
        slots = self.slots(k)[0]
        dirs = self.windings[k]['dir']
        # turns = self.windings[k]['N']
        r = len(slots)//len(dirs)
        curr = np.concatenate([np.array(dirs)*(1 - 2*(n % 2))
                               for n in range(r)])

        NY = 4096
        y = np.zeros(NY*self.Q//t)
        for i in range(1, self.Q//t+1):
            if i in set(slots):
                y[NY*(i-1)+NY//2] = np.sum(curr[slots == i])
        yy = [np.sum(y[:i+1]) for i in range(0, len(y))]
        yy[:NY//2] = yy[-NY//2:]
        yy = np.tile(yy-np.mean(yy), t)
        yy /= np.max(yy)
        # y = np.tile(y,t)

        N = len(yy)
        Y = np.fft.fft(yy)
        imax = np.argmax(np.abs(Y[:N//2]))
        a = 2*np.abs(Y[imax])/N
        freq = np.fft.fftfreq(N, d=taus/NY)
        T0 = np.abs(1/freq[imax])
        alfa0 = np.angle(Y[imax])
        # if alfa0 < 0: alfa0 += 2*np.pi
        pos_fft = np.linspace(0, self.Q/t*taus, self.p//t*60)
        D = (a*np.cos(2*np.pi*pos_fft/T0+alfa0))
        return dict(
            pos=[i*taus/NY for i in range(len(y))],
            mmf=yy[:NY*self.Q//t].tolist(),
            alfa0=-alfa0/self.p,
            pos_fft=pos_fft.tolist(),
            nue=np.arange(0, 9*self.p).tolist(),
            mmf_nue=(2*np.abs(Y[:9*self.p])/N).tolist(),
            mmf_fft=D.tolist())

    def zoneplan(self):
        taus = 360/self.Q
        dphi = 1e-3
        slots = {k: [s-1 for s in self.slots(k)[0]]
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

        z = ([[d*s for s, d in zip(u, ud)] for u, ud in zip(upper, udirs)],
             [[d*s for s, d in zip(l, ld)] for l, ld in zip(lower, ldirs)])
        # complete if not  basic winding:
        Qb = self.Q//num_basic_windings(self.Q, self.p, self.l)
        if max([abs(n) for m in z[0] for n in m]) < Qb:
            return [[k + [-n+Qb//2 if n < 0 else -(n+Qb//2) for n in k]
                     for k in m] for m in z]
        return z

    def diagram(self) -> ET.Element:
        """return winding diagram as svg element"""
        coil_len = 25
        coil_height = 4
        dslot = 8
        arrow_head_length = 2
        arrow_head_width = 2
        strokewidth = [f"{w}px" for w in [0.25, 0.5]]

        z = self.zoneplan()
        xoff = 0
        if z[-1]:
            xoff = 0.75
        yd = dslot*self.yd
        mh = 2*coil_height/yd
        slots = sorted([abs(n) for m in z[0] for n in m])
        smax = slots[-1]*dslot
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        svg = ET.Element("svg", dict(
            version="1.1", xmlns="http://www.w3.org/2000/svg",
            viewBox=f"0, -30, {slots[-1] * dslot + 15}, 40"))
        g = ET.SubElement(svg, "g", {"id": "teeth", "fill": "lightblue"})
        for n in slots:
            e = ET.SubElement(g, "rect", {
                "x": f"{n * dslot + dslot/4}",
                "y": f"{-coil_len + 1}",
                "width": f"{dslot/2}",
                "height": f"{coil_len - 2}"})

        g = ET.SubElement(svg, "g",
                          {"id": "labels",
                           "text-anchor": "middle",
                           "dominant-baseline": "middle",
                           "style": "font-size: 0.15em; font-family: sans-serif;"})
        for n in slots:
            t = ET.SubElement(g, "text", {
                "x": f"{n*dslot}",
                "y": f"{-coil_len / 2}"}).text = str(n)

        g = ET.SubElement(svg, "g", {"id": "coils",
                                     "fill": "none",
                                     "stroke-linejoin": "round",
                                     "stroke-linecap": "round"})

        for i, layer in enumerate(z):
            b = -xoff if i else xoff
            w = i if self.yd > 1 else 0
            for m, mslots in enumerate(layer):
                for k in mslots:
                    slotpos = abs(k) * dslot + b
                    pc = [f"L {slotpos} {-coil_len//2+2} M {slotpos} {-coil_len//2-1}",
                          f"L {slotpos} {-coil_len}"]
                    if (i == 0 and (k > 0 or (k < 0 and self.l > 1))):
                        # first layer, positive dir or neg. dir and 2-layers:
                        #   from right bottom
                        if slotpos + yd > smax+b:
                            dx = dslot if yd > dslot else yd/4
                            ph = [f"M {slotpos+yd//2-xoff+dx} {coil_height-mh*dx}",
                                  f"L {slotpos+yd//2-xoff} {coil_height} L {slotpos} 0"]
                            pt = [f"L {slotpos+yd//2-xoff} {-coil_len-coil_height}",
                                  f"L {slotpos+yd//2-xoff+dx} {-coil_len-coil_height+mh*dx}"]
                        else:
                            ph = [
                                f"M {slotpos+yd//2-xoff} {coil_height} L {slotpos} 0"]
                            pt = [
                                f"L {slotpos+yd//2-xoff} {-coil_len-coil_height}"]
                    else:
                        # from left bottom
                        if slotpos - yd < 0:  # and slotpos - yd > -3*dslot:
                            dx = dslot if yd > dslot else yd/4
                            ph = [f"M {slotpos-yd//2+xoff-dx} {coil_height-mh*dx}",
                                  f"L {slotpos-yd//2+xoff} {coil_height} L {slotpos} 0"]
                            pt = [f"L {slotpos-yd//2+xoff} {-coil_len-coil_height}",
                                  f"L {slotpos-yd//2+xoff-dx} {-coil_len-coil_height+mh*dx}"]
                        else:
                            ph = [
                                f"M {slotpos-yd//2+xoff} {coil_height} L {slotpos} 0"]
                            pt = [
                                f"L {slotpos-yd//2+xoff} {-coil_len-coil_height}"]
                    e = ET.SubElement(g, "path", {
                        "d": ' '.join(ph + pc + pt),
                        "stroke-width": strokewidth[w],
                        "stroke": coil_color[m]})

        for i, layer in enumerate(z):
            for m, mslots in enumerate(layer):
                for k in mslots:
                    x = abs(k) * dslot
                    if i:
                        x -= xoff
                    else:
                        x += xoff
                    if k > 0:
                        y = coil_len * .88
                        points = [
                            (x, -y),
                            (x - arrow_head_width / 2, -y + arrow_head_length),
                            (x + arrow_head_width / 2, -y + arrow_head_length)]
                    else:
                        y = coil_len * .12
                        points = [
                            (x, -y),
                            (x - arrow_head_width / 2, -y - arrow_head_length),
                            (x + arrow_head_width / 2, -y - arrow_head_length)]
                    ET.SubElement(svg, "polygon", {
                        "points": " ".join([f"{x},{y}" for (x, y) in points]),
                        "fill": f"{coil_color[m]}",
                        "stroke": "none"})

        return ET.tostring(svg, encoding='unicode')

    def write(self, name, workdir='.'):
        """creates WID file"""
        import pathlib
        with open(pathlib.Path(workdir) / (name + '.WID'), 'w') as fp:
            if 'slots' in self.windings[1]:
                inp = sorted([(k if s > 0 else -k, abs(s), l, n)
                              for k in self.windings
                              for s, l, n in zip(
                    self.windings[k]['slots'],
                    self.windings[k]['layer'],
                    self.windings[k]['N'])],
                    key=lambda x: (x[1], x[2]))
                fp.write('\n'.join([
                    f'Windings input data: {name}',
                    ' Number of coil sides:',
                    f'          {len(inp)}',
                    'Type of machine: 1 = Rot, 21 = Lin-x, 22 = Lin-y',
                    '           1',
                    'Index  w-keys     N-turns         Layer', '']))

                for i, t in enumerate(inp):
                    fp.write(f'{i+1:5d} {t[0]:8d} {t[3]:10.3f} {t[2]:10d}\n')
            else:
                inp = sorted([(k if d > 0 else -k, r, phi, n)
                              for k in self.windings
                              for d, r, phi, n in zip(
                    self.windings[k]['dir'],
                    self.windings[k]['R'],
                    self.windings[k]['PHI'],
                    self.windings[k]['N'])],
                    key=lambda x: (x[1], x[2]))
                fp.write('\n'.join([
                    f'Windings input data: {name}',
                    ' Number of coil sides:',
                    f'          {len(inp)}',
                    'Type of machine: 1 = Rot, 21 = Lin-x, 22 = Lin-y',
                    '           1',
                    'Index  w-keys     N-turns         R[mm]     PHI[deg]', '']))
                for i, t in enumerate(inp):
                    fp.write(
                        f'{i+1:5d} {t[0]:8d} {t[3]:10.3f} {t[1]:8.3f} {t[2]:8.3f}\n')
            fp.write('\n'.join([
                'Number of windings saved :',
                f'         {len(self.windings)}',
                'W-Key Coil-Current [A]   W-Types: (=1 :wire&cur)'
            ] + [
                f'   {k}          0.00000000       0.00000000                            1'
                for k in self.windings
            ] + ['   0', '']))


if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    from xml.etree import ElementTree as ET

    if sys.argv[1:]:
        bch = femagtools.bch.read(sys.argv[1])
        wdgs = Winding(bch)
    else:
        testdata = [
            dict(Q=90, p=12, m=3, l=2, coilwidth=1),
            dict(Q=54, p=6, m=3, l=2, coilwidth=5),
            dict(Q=168, p=7, m=3, l=2, coilwidth=10)]

        wdgs = Winding(testdata[1])

    c = wdgs.mmf()
    # print('alfa0={0:6.3f}'.format(wdgs.axis()/np.pi*180))

    plt.title('Q={0}, p={1}, alfa0={2:6.3f}'.format(
        wdgs.Q, wdgs.p, c['alfa0']/np.pi*180))
    plt.plot(np.array(c['pos'])/np.pi*180, c['mmf'])
    plt.plot(np.array(c['pos_fft'])/np.pi*180, c['mmf_fft'])

    phi = [c['alfa0']/np.pi*180, c['alfa0']/np.pi*180]
    y = [min(c['mmf_fft']), 1.1*max(c['mmf_fft'])]
    plt.plot(phi, y, '--')
    plt.annotate("", xy=(phi[0], y[0]),
                 xytext=(0, y[0]), arrowprops=dict(arrowstyle="->"))

    plt.grid()
    plt.show()

    svg = wdgs.diagram()
    ET.ElementTree(svg).write('wind.svg')
    print('SVG file "wind.svg" created')
