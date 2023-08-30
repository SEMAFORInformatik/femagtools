"""
    femagtools.ntib
    ~~~~~~~~~~~~~~~

    NTIB / LOS files handling



"""
import logging
import math

logger = logging.getLogger(__name__)


def create(speed, current, beta, r1=0, m=3):
    """return Ntib info"""
    return ['Speed[1/min], Current[A], beta[Degr], R_STATOR[Ohm], n-p',
            'Number of values:      {0}'.format(
                min(len(speed), len(current), len(beta)))] + \
                ['{0:12.1f}{1:12.3f}{2:12.3f}{3:12f}{4:12f}'.format(
                    60*n, math.sqrt(2.0)*i1,
                    b,
                    r1, m)
                 for n, i1, b in zip(speed, current, beta)]


def toFloat(s, fac=1.0):
    try:
        return float(s)*fac
    except ValueError:
        return float('nan')




def read_los_content(content):
    """return dict of losses in LOS-file content"""
    result = dict(speed=[],
                  torque=[],
                  i1=[],
                  beta=[],
                  stajo=[],
                  staza=[],
                  rotfe=[],
                  magnet=[],
                  winding=[],
                  total=[])
    started = False
    for l in content:
        if not started and l.startswith('[1/min]   '):
            started = True
        elif started:
            r = l.strip().split()
            if len(r) > 7:
                result['speed'].append(toFloat(r[0], 1./60))
                result['torque'].append(toFloat(r[1]))
                result['i1'].append(toFloat(r[2]))
                result['beta'].append(toFloat(r[3]))
                pfe1 = toFloat(r[4])
                pfe2 = 0
                p = 4
                if len(r) > 8:
                    pfe2 = toFloat(r[5])
                    p = 5
                result['stajo'].append(pfe1)
                result['staza'].append(pfe2)
                result['rotfe'].append(toFloat(r[p+1]))
                result['magnet'].append(toFloat(r[p+2]))
                result['winding'].append(toFloat(r[p+3]))

                try:
                    result['stafe'] = result['stajo'] + result['staza']
                    result['total'].append(sum([result[k][-1]
                                                for k in ('stajo',
                                                          'staza',
                                                          'rotfe',
                                                          'magnet',
                                                          'winding')]))
                except KeyError:
                    result['total'].append(None)
    logger.info("num rows %d", len(result['total']))
    return result


def read_los(filename):
    """return dict of losses in LOS-file"""
    logger.info("read loss file: %s", filename)
    with open(filename) as f:
        return read_los_content(f.readlines())

    # empty
    return dict(speed=[],
                torque=[],
                i1=[],
                beta=[],
                stajo=[],
                staza=[],
                rotfe=[],
                magnet=[],
                winding=[],
                total=[])
