import numpy as np

def fft(pos, y, pmod=0):
    """calculate fft spectrum of y and return samples,
    values, amplitude and phase of base harmonic

    Arguments:
      pos: (list of floats) sample positions
      y: (list of floats) y values
      pmod: number of poles in model (ignored if 0)
    """
    model_angle = pos[-1] - pos[0]
    ntiles = int(round(360/model_angle))

    if pmod:
        negative_periodic = pmod % 2
    else:
        negative_periodic = np.abs(y[0] - y[-1])/np.max(y) > 1

    if negative_periodic:
        yx = np.append(
            np.concatenate(
                [n*y[:-1]
                 for n in [m % 2 or -1
                           for m in range(1, ntiles+1)]]),
            y[0])
    else:
        yx = np.append(
            np.tile(y[:-1], ntiles),
            y[0])

    N = len(yx)
    # compute DFT from induction (eliminate DC offset)
    a0 = np.mean(yx)
    Y = np.fft.fft(yx-a0)

    # find the peak (amplitude of base harmonic)
    i = np.argmax(np.abs(Y[:N//2]))

    a = 2*np.abs(Y[i])/N
    freq = np.fft.fftfreq(N, d=pos[1]-pos[0])
    nmax = min(18*ntiles, N//2)
    T0 = 0
    if abs(freq[i]) > 0:
        T0 = np.abs(1/freq[i])
        npoles = 2*int(np.ceil(360/T0))
        nmax = min(9*npoles, N//2)

    alfa0 = np.angle(Y[i])

    return {'a': a, 'a0': a0, 'T0': T0, 'alfa0': alfa0,
            'nue': (2*np.abs(Y[:nmax])/N).tolist(),
            'yi': yx.tolist()}
