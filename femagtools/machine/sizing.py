"""
  femagtools.machine.sizing
  ~~~~~~~~~~~~~~~~~~~~~~~~~

  general design of an AC electrical machine

"""
import numpy as np
import logging

logger = logging.getLogger("femagools.machineo.sizing")


PM_DEFAULTS = dict(
    airgap=1.5e-3,  # airgap width m
    eta=0.92,  # efficiency
    cos_phi=0.95,  # power factor
    m=3,  # number of phases
    ui_u=0.62,  # U ind / U
    lda=1.33,  # length/taup ratio 0.5 .. sqrt(p)
    J=3.8e6,  # current density 3 .. 6 A/mm²
    sigmas=17e3,  # shear force 10 .. 45 kN/m2
    Ba=0.7,  # flux density in airgap 0.6 .. 1.5 T
    Bth=1.5,  # flux density in teeth 1.5 .. 2 T
    By=1.2,  # flux density in yoke 1.2 .. 1.5 T
    kq=0.42,  # stator winding fill factor 0.35 .. 0.75
    mag_width=0.8,  # rel magnet width 0.6 .. 1
    Hc=700,  # max. coercitive field strength, 500 .. 900 kA/m
    brem=1.15,  # remanence 0.3 .. 1.3 T
    demag=1.5,  # safety factor for demagnetisation
    external_rotor=False,
    coil_span=0,
    hfe=1e-3  # iron height between magnet and airgap (IPM)
)

SM_DEFAULTS = dict(
    airgap=1.5e-3,  # airgap width m
    eta=0.92,  # efficiency
    cos_phi=0.95,  # power factor
    m=3,  # number of phases
    ui_u=0.62,  # U ind / U
    lda=1.33,  # length/taup ratio 0.5 .. sqrt(p)
    J=3.8e6,  # current density 3 .. 6 A/mm²
    sigmas=17e3,  # shear force 10 .. 45 kN/m2
    Ba=0.7,  # flux density in airgap 0.6 .. 1.5 T
    Bth=1.5,  # flux density in stator teeth 1.5 .. 2 T
    By=1.2,  # flux density in stator yoke 1.2 .. 1.5 T
    Bthr=1.7,  # flux density in rotor teeth 1.5 .. 2 T
    Byr=1.4,  # flux density in rotor yoke 1.2 .. 1.5 T
    kq=0.42,  # stator winding fill factor 0.35 .. 0.75
    kqr=0.6,  # rotor winding fill factor 0.45 .. 0.7
    external_rotor=False,
    coil_span=0
)

IM_DEFAULTS = dict(
    airgap=0.75e-3,  # airgap width m
    eta=0.87,  # efficiency
    cos_phi=0.8,  # power factor
    m=3,  # number of phases
    ui_u=0.62,  # U ind / U
    lda=1.33,  # length/taup ratio 0.5 .. sqrt(p)
    J=3e6,  # current density 3 .. 6 A/mm²
    sigmas=14e3,  # shear force 10 .. 45 kN/m2
    Ba=0.7,  # flux density in airgap 0.6 .. 1.5 T
    Bth=1.5,  # flux density in teeth 1.5 .. 2 T
    By=1.2,  # flux density in yoke 1.2 .. 1.5 T
    kq=0.42,  # winding fill factor 0.35 .. 0.75
    kqr=.9,   # bar fill factor
    external_rotor=False,
    coil_span=0
)


class ValidationError(Exception):
    pass


def iron_losses(dy, lfe, hy, wt, ht, q1, B, f):
    dens = 7.8e3  # mass density of iron kg/m³
    ch = 3
    ce = 2
    Bo = 1.5
    fo = 50
    Wfe = ch*(B/Bo)**1.8*f/fo + ce*(B/Bo)**2*(f/fo)**2
    vol = lfe*(np.pi*(dy**2/4 - (dy/2-hy)**2) + q1*(ht*wt))
    return Wfe * dens*vol


def wdg_resistance(w1, l, d):
    S = 56e6  # conductivity of copper 1/Ohm m
    a = np.pi*d**2/4
    return 2*w1*l/S/a


def current_losses(w1, l, d, i1):
    # w1: num windings/phase
    # l: wire length/slot
    # d: wire diameter
    return wdg_resistance(w1, l, d) * i1**2


def check_symmetry_conditions(Q1, p, layers, m):
    # if Q1 % 2*p == 0:
    #    raise ValidationError("Number of slots Q and poles 2*p are not prime to each other")

    K = 2 if layers == 1 else 1
    urwick = np.gcd(Q1//K, p)
    if Q1 % m != 0:
        raise ValidationError(
            "Number of slots Q and phases are not consistent")

    K = 2 if m % 2 == 0 else 1
    if Q1 % K*m*urwick != 0:
        raise ValidationError(
            "Number of slots Q and phases are not consistent")


# number of rotor slots of IM:
# recommendations (Alger, Kuhlmann, Veinott):
# p, Q1 --> Qr
Q2R = {
    1: {36: [26, 28, 44],
        48: [38, 40, 56],
        54: [54],
        60: [52, 68, 78]},
    2: {36: [26, 28, 44],
        48: [34, 38, 40, 56],
        60: [34, 44, 46, 76],
        72: [58]},
    3: {36: [46, 48],
        48: [58, 60, 64, 68],
        54: [42, 66],
        72: [54, 58, 84, 88]},
    4: {36: [48, 52],
        48: [58, 64],
        54: [70],
        72: [58, 88]}}
# Richter:
Q2RR = {
    1: {24: [28, 16, 22],
        36: [24, 28, 48, 16],
        48: [40, 52],
        60: [48]},
    2: {36: [24, 40, 42, 60, 30, 44],
        48: [60, 84, 56, 44],
        60: [72, 48, 84, 44]},
    3: {36: [42, 48, 54, 30],
        54: [72, 88, 48],
        72: [96, 90, 84, 54]},
    4: {36: [48],
        48: [72, 60],
        72: [96, 84]}}

Q2IEC = {
    1: {18: [22],
        24: [18, 22],
        36: [],
        48: [],
        54: [],
        60: [],
        72: []}
}


def _im_stator_slots(p):
    return [k for k in Q2R[p].keys()]


def _rotor_slots(Q1, p):
    """returns a list of valid rotor slot numbers"""
    # Slots Qr of Induction Machines:
    #
    #  Qr <= 1.25 Qs
    #  Qr != Qs +/- 2p
    #  Qr != Qs
    #  Qr != Qs +/- 1
    # Slot area Acr: 0.7 - 0.9 Acs
    #
    # Examples of 4 poles SCIM Qs, Qr:
    # Audi e-Tron 48 58
    # Tesla S, X  60 74
    # Renault Twizy 36 48
    # others: 24/22, 36/24
    try:
        return Q2R[p][Q1]
    except KeyError:
        pass
    d = 0.25
    sync = Q1+p, Q1-p, Q1+2*p, Q1-2*p, Q1+5*p, Q1-5*p
    noise = Q1+1, Q1-1, Q1+2, Q1 + p+1, Q1-p-1, Q1+p+2, Q1-p-2
    return sorted(set(list(range(int((1-d)*Q1), int((1+d)*Q1)))) -
                  set((Q1, ) + sync + noise))


def _stator_slots(par, slots):
    """find a suitable number of slots
    (condition: 3 < slot height / width  < 5)"""
    q = []
    p = par['p']
    m = 3
    if slots:
        slotset = slots
    else:
        slotset = [n for n in range(12, 110, 2*m)
                   if n % (m*np.gcd(n, p)) == 0]
    for Q1 in slotset:
        par['Q1'] = Q1
        r = get_stator_dimensions(par)
        hb = r['hns']/r['bns']
        if hb > 5:
            break
        if hb > 3:
            q.append((Q1, Q1//np.gcd(Q1, p), abs(4-hb), r['kw']))

    if q:
        qhbmin = np.argmin(q, axis=0)
        logger.debug("q %s", q)
        # check sim factor, height/width ratio, winding factor
        if qhbmin[1] == len(q)-1:
            return q[qhbmin[1]][0]
        elif q[qhbmin[1]][2] < q[qhbmin[1]+1][2]:
            return q[qhbmin[1]][0]
        return q[qhbmin[1]+1][0]
    return slotset[0]


def get_stator_dimensions(par, slots=[]):
    # nominal (rated) operating parameters
    pnom = par['pnom']
    speednom = par['speed']
    tnom = pnom/(2*np.pi*speednom)
    p = par['p']  # pole pairs
    lda = par['lda']  # length/taup ratio 0.5 * np.sqrt(p)
    sigmas = par['sigmas']  # shear force, typically 10 .. 45 kNm/m²

    # airgap diameter
    Da = (4*p * tnom/(np.pi**2 * lda * sigmas)) ** (1./3)

    # pole width
    taup = np.pi * Da/(2*p)
    lfe = taup*lda

    if 'udc' in par:
        u1nom = 0.9*par['udc']/np.sqrt(2)/np.sqrt(3)  # phase voltage
    else:
        u1nom = par['u1']  # phase voltage

    eta = par['eta']  # efficiency
    cos_phi = par['cos_phi']

    # winding and poles
    m = par['m']  # num phases
    if 'Q1' not in par:
        Q1 = _stator_slots(par, slots)
    else:
        Q1 = par['Q1']
    coil_span = par.get('coil_span', 0)
    if np.gcd(Q1, 2*p) == coil_span:
        coil_span = 0
    layers = 1 if coil_span == 0 else 2
    req_poles = np.gcd(Q1, 2*p)
    if req_poles != 2*p:
        layers = 2

    check_symmetry_conditions(Q1, p, layers, m)

    Ui_U = par['ui_u']

    # design parameters
    # AJ = papyr['AJ'] # 100*1e9 # thermal load, typically 100 .. 300 A²/mm³
    J = par['J']   # current density 3 .. 6 A/mm²
    Ba = par['Ba']     # airgap flux density, typically 0.6 .. 1.5 T
    kq = par['kq']      # winding fill factor, typically 0.35 .. 0.75
    # flux density in teeth, typically 1.7 .. 2.2 T
    Bth = par['Bth']
    By = par['By']     # flux density in yoke, typically 1.2 .. 1.5 T

    # winding factor
    Nt = Q1/np.gcd(Q1, p)
    qt = Nt / (2*m)  # if Nt % 2 == 0 else Nt / m
    kwz = np.sin(np.pi/(2*m))/(qt*np.sin(np.pi/(2*m*qt)))
    kws = np.sin(np.pi*p*coil_span/Q1) if coil_span != 0 else 1.0
    kw = kwz * kws

    # flux density amplitude of base harmonic
    if 'brem' in par:
        mag_width = par['mag_width']
        Bd1 = 4.0/np.pi*Ba*np.sin(np.pi/2.0*mag_width)
    else:
        Bd1 = Ba
   # ---------------------------------------------
    f1nom = speednom*p

    # stator windings
    psi1 = 2.0/np.pi*taup*lfe*Bd1

    Ui = Ui_U * u1nom
    N = round(np.sqrt(2)*Ui/(2*np.pi*f1nom*kw*psi1), 1)

    # current loading
    # A = np.sqrt(2)*sigmas/kw/Ba
    I1 = pnom/(m*eta*cos_phi*u1nom)
    A = 2*m*N*I1/np.pi/Da

    # slot area
    # J = AJ/A
    hs1 = 1e-3  # slot opening height
    taus = np.pi/Q1
    ans = taus*Da*A/(kq*J)
    bds = taus*(Da+2*hs1)*Bd1/Bth
    bns = taus*(Da+2*hs1) - bds

    hns = (-bns + np.sqrt(bns**2 + 4*ans*np.tan(taus)))/2/np.tan(taus)
    hys = psi1/2/lfe*By

    # airgap and yoke diameter
    airgap = par['airgap']
    if par['external_rotor']:
        Da1 = Da - 2*airgap
        Dy1 = Da1 - 2*(hys + hns + hs1)
        Da2 = Da
    else:
        Da1 = Da
        Dy1 = Da1 + 2*(hys + hns + hs1)
        Da2 = Da1-2*airgap
    r = dict(
        lfe=round(lfe, 3),
        Dy1=round(Dy1, 3),
        Da1=round(Da1, 4),
        Da2=round(Da2, 4),
        ans=round(ans, 6),
        hns=round(hns, 4),
        bns=round(bns, 4),
        A=round(A, 3),
        AJ=round(J*A, 0),
        ess=round(1e-3*pnom/(60*speednom)/(Da1**2*lfe), 4),
        trv=round(1e-3*pnom/(2*np.pi*speednom)/(np.pi*Da1**2/4*lfe), 1),
        w1=int(N),
        kw=round(kw, 4),
        q=qt,
        i1=round(I1, 3),  # np.pi*Da1*A/2/m/N
        psi1=round(psi1, 5))

    slotwidth = np.floor(0.3*np.pi*Da1/Q1*2000+0.5)/2000.0
    tw = bds  # np.pi*(Da1+hns)/Q1 - bns
    toothwidth = np.floor(tw*2000+0.5)/2000.0
    bns1 = np.floor((np.pi*Da1/Q1 - tw)*2000.0+0.5)/2000.0
    bns2 = np.floor(((np.pi*Da1+2*hns)/Q1 - tw)*2000+0.5)/2000.0

    r['pfe'] = round(iron_losses(Dy1, lfe, (Dy1/2-Da1/2-hns),
                                 toothwidth, hns, Q1,
                                 Bd1, f1nom), 1)

    if coil_span == 0:
        tp = Q1/2./p
        q = tp/m
        coil_span = Q1/p/2 if q >= 1 else 1

    middle_line = 1 if coil_span == 1 else 2

    # stator 3
    r['stator'] = dict(
        u1nom=round(u1nom, 1), f1nom=round(f1nom, 1),
        num_slots=Q1,
        nodedist=1,
        # num_slots_gen = req_poles*Q1/2/p,
        statorRotor3=dict(
            slot_width=slotwidth,
            tooth_width=toothwidth,
            slot_height=hs1+r['hns'],
            slot_top_sh=1.0,
            slot_h1=hs1,
            slot_h2=hs1,
            slot_r1=1.5e-3,
            slot_r2=1.5e-3,  # bns2/2,
            wedge_width1=0,  # bns1,
            wedge_width2=0,
            middle_line=0 if layers < 2 else middle_line))
    # coils per phase and pole
    # q = Q/2/p/m
    # num wires per coil side
    # n = a*N / 2 / p / q
    num_wires = round(2*m*N/layers/Q1)
    dwire = 2*np.sqrt(ans*kq/layers/np.pi/num_wires)
    relculen = 1.4
    r1 = wdg_resistance(N, relculen*lfe, dwire)
    r['pcu'] = round(m*r1*I1**2, 1)

    r['windings'] = dict(
        wire_diam=round(dwire, 5),
        num_phases=m,
        cufilfact=kq,
        culength=1.4,
        num_layers=layers,
        resistance=round(r1, 4),
        coil_span=int(coil_span),
        num_wires=int(num_wires))

    return r


def _get_magnet_height(I1, N, kw, par):
    airgap = par['airgap']
    Ba = par['Ba']
    p = par['p']
    m = par['m']

    Hc = par['Hc']
    Hc = Hc*1e3  # unit kA/m
    # Safety Factor for demagnetization
    demag = par['demag']
    I1max = 4*I1
    THETA1 = m/np.pi*np.sqrt(2)*I1max*N*kw/p
    # minimal magnet height with safety factor
    hM = THETA1/Hc*demag
    Br = par['brem']
    if hM/(hM+airgap)*Br < Ba:
        raise ValidationError(
            f'Flux density in airgap {Ba:4.2f} T is too high, must be smaller than {hM/(hM+airgap)*Br:4.2f} T')
    return hM


def _get_magnet_yoke_diameter(psi1, lfe, Da2, hM, By, external_rotor):
    hys = psi1/(2*By*lfe)
    if external_rotor:
        return (Da2 + 1.2*2*hM + 2*hys)

    return (Da2 - 1.2*2*hM - 2*hys)


def get_surface_magnet_dimensions(I1, N, kw, psi1, lfe, Da2, par):
    hM = _get_magnet_height(I1, N, kw, par)
    Dy2 = _get_magnet_yoke_diameter(
        psi1, lfe, Da2, hM, par['By'], par['external_rotor'])
    return dict(
        Da2=np.floor(Da2*1e4+0.5)/1e4,
        Dy2=np.floor(Dy2*1e4+0.5)/1e4,
        hM=round(hM, 5),
        nodedist=1.0,
        magnetSector=dict(
            magn_height=np.floor(hM*2000+0.5)/2000.0,
            magn_width_pct=par['mag_width'],
            condshaft_r=round(Dy2/2, 4),
            magn_num=1,
            magn_rfe=0.0,
            magn_len=1.0,
            magn_shape=0.0,
            bridge_height=0.0,
            bridge_width=0.0,
            magn_ori=2,
            magn_type=1))


def get_interior_magnet_dimensions(I1, N, kw, psi1, lfe, Da2, par):
    r = dict()
    hM = round(_get_magnet_height(I1, N, kw, par), 5)
    Dy2 = _get_magnet_yoke_diameter(
        psi1, lfe, Da2, hM, par['By'], par['external_rotor'])

    kr = 0.9
    # pole pitch angle
    alphap = np.pi/par['p']
    hfe = par['hfe']
    tap = np.tan(alphap/2)
    sip = np.sin(alphap/2)
    a = (1+1/tap**2)
    b = 4/tap*(hfe/sip+hM)
    c = -Da2**2 + 4*Da2*hfe - 4*hfe**2 + 4*(hfe/sip+hM)**2
    wM = (-b+np.sqrt(b**2 - 4*a*c))/2/a
    hp = Da2/2 - (wM/tap/2 + hfe/sip + hM)

    # Magnet iron 4 BFE
    mag_dm_min = Da2/2 - np.sqrt((Da2**2 - wM**2)/4)
    mag_dm_max = (Da2-Dy2)/2 - hM
    magn_di_ra = mag_dm_min * kr + mag_dm_max * (1 - kr)
    mag_corner_center_distance = np.sqrt(
        (Da2/2-magn_di_ra-hM)**2 + (wM/2)**2)
    iron_bfe = (mag_corner_center_distance *
                np.sin(alphap/2 -
                       np.arcsin(wM/2/mag_corner_center_distance)))
    return dict(
        Da2=np.floor(Da2*1e4+0.5)/1e4,
        Dy2=np.floor(Dy2*1e4+0.5)/1e4,
        hM=hM,
        nodedist=1.0,
        magnetIron=dict(
            magn_height=np.floor(hM*2000+0.5)/2000.0,
            magn_width=round(wM, 4),
            air_triangle=1,
            magn_rem=par['brem'],
            iron_height=round(hfe, 4),
            gap_ma_iron=0,
            bridge_height=0.0,
            bridge_width=0.0,
            magn_ori=2,
            condshaft_r=round(Dy2/2, 4),
            iron_shape=Da2/2))


def get_im_rotor_dimensions(A, Da2, psi1, lfe, par, rtype='rotorKs2'):
    r = dict()
    r['Da2'] = Da2
    if 'Q2' not in par:
        r['num_slots'] = _rotor_slots(par['Q1'], par['p'])[0]
    else:
        r['num_slots'] = par['Q2']
    # slot area
    # J = AJ/A
    J = par['J']
    hs1 = 1e-3
    alfar = np.pi/r['num_slots']
    kq = par['kqr']
    Ar = alfar*Da2*A/J/kq
    wt = alfar*(Da2-2*hs1)*par['Ba']/par['Bth']
    wr = alfar*(Da2-2*hs1) - wt
    #hr = (wr + np.sqrt(wr**2 - 4 * Ar*np.tan(alphar)))(2*np.tan(alphar))
    hr = (wr + np.sqrt(wr**2 +
                       4*Ar*(1-2*np.tan(alfar))))/2/(1-2*np.tan(alfar))
    taup = np.pi * Da2/(2*par['p'])
    hyr = psi1/2/lfe*par['By']
    r['Dy2'] = round(Da2 - 2*hr - 2*hyr, 4)
    logger.info("Dy2 %f Da2 %f hys %f hr %f",
                r['Dy2']*1e3, Da2*1e3, hyr*1e3, hr*1e3)
    slotwidth = 1e-3
    Q2 = r['num_slots']
    r1 = wr/2-slotwidth
    r2 = (Da2/2-hr-hs1)*np.tan(alfar)
    if rtype == 'statorRotor3':
        r['statorRotor3'] = dict(
            slot_width=slotwidth,
            tooth_width=round(wt, 4),
            slot_height=round(hr+r2, 4),
            slot_top_sh=1.0,
            slot_h1=round(hs1, 4),
            slot_h2=round(hs1+r1, 4),
            slot_r1=round(r1, 4),
            slot_r2=round(r2),
            wedge_width1=0,
            wedge_width2=0,
            middle_line=0)
    elif rtype == 'rotorAsyn':
        r['rotorAsyn'] = dict(
            slot_bs2=0.1e-3,
            slot_hs2=0.5e-3,
            slot_b32=0.0,
            slot_h32=0.0,
            slot_b42=0.0,
            slot_h42=0.0,
            slot_b52=round(wr, 4),
            slot_b62=3e-3,
            slot_h52=2.5e-3,
            slot_h62=round(hr, 4),
            slot_h72=2e-3)
    else:
        r['rotorKs2'] = dict(
            slot_angle=round(2*alfar*180/np.pi, 2),
            slot_height=round(hr+r1+r2, 4),
            slot_topwidth=round(wr, 4),
            slot_width=slotwidth,
            slot_h1=hs1,
            slot_h2=0,
            slot_r1=1e-3,  # r1,
            slot_r2=1e-3,  # r2,
            middle_line=0)

    return r


def get_sm_rotor_dimensions(A, psi1, lfe, Da2, par):
    r = dict()
    r['Da2'] = Da2
    r['num_wires'] = 1

    alphap = np.pi/par['p']
    wp = 0.9*Da2*np.sin(alphap/2)
    hp = Da2/2*(1-np.cos(alphap/2))
    wc = alphap/2*Da2*par['Ba']/par['Bthr']
    hfmax = (wp-wc)/2/np.tan(alphap/2)
    mue0 = 4*np.pi*1e-7

    wf = wp - wc
    anr = par['airgap']*par['Ba']/mue0/par['J']/par['kqr']

    a = np.tan(alphap/2)/2
    b = (wp/2 - wc + np.tan(alphap/2)*(Da2/2 - hp))/2
    hc = min(max(4*anr/(wc+wp), 0.75*hfmax), hfmax)

    hyr = psi1/2/lfe*par['Byr']
    pr = 0.95*Da2/2

    r['Dy2'] = round(Da2 - 2*hp - 2*hyr - 2*hc, 4)
    if r['Dy2'] < 6e-3:
        hyr = hyr + r['Dy2']/2
        r['Dy2'] = 6e-3  # minimal size

    r['rot_hsm'] = dict(
        gap_pol_shaft=0.0,
        core_height=round(hc, 4),
        pole_height=round(hp, 4),
        pole_rad=round(pr, 4),
        core_width2=round(wc, 4),
        core_width1=round(wc, 4),
        pole_width_r=round(wp, 4),
        pole_width=round(wp, 4),
        slot_width=0,
        slot_height=0,
        damper_diam=0,
        damper_div=0
    )
    return r


def _set_defaults(par, defaults):
    for k in defaults:
        if k not in par:
            par[k] = defaults[k]


def _set_pm_defaults(par):
    _set_defaults(par, PM_DEFAULTS)


def _set_im_defaults(par):
    _set_defaults(par, IM_DEFAULTS)


def _set_sm_defaults(par):
    _set_defaults(par, SM_DEFAULTS)


def _set_genpars(r, poles):
    r.pop('Da2')
    for k in ('magnet', 'rotor'):
        if k in r:
            break
    r['outer_diam'] = r.pop('Dy1')
    r['inner_diam'] = r[k].pop('Dy2')
    r['airgap'] = (r['Da1'] - r[k].pop('Da2'))/2
    r['bore_diam'] = r.pop('Da1')
    r['poles'] = poles
    for k in 'u1nom', 'f1nom':
        r[k] = r['stator'].pop(k)


def spm(pnom, speed, p, **kwargs):
    """returns dimension of a SPM machine"""
    par = dict(
        pnom=pnom, speed=speed, p=p)
    par.update(kwargs)

    _set_pm_defaults(par)
    # stator and magnet parameters of surface mounted magnet machine
    r = get_stator_dimensions(par)

    # magnet parameters
    r['magnet'] = get_surface_magnet_dimensions(
        r['i1'], r['w1'], r['kw'], r['psi1'], r['lfe'], r['Da2'],
        par)

    r['magnet'].pop('hM')
    _set_genpars(r, 2*par['p'])
    r['name'] = f"SPM-{r['poles']}"

    return r


def ipm(pnom, speed, p, **kwargs):
    """returns dimension of a IPM machine"""
    par = dict(
        pnom=pnom, speed=speed, p=p)
    par.update(kwargs)
    _set_pm_defaults(par)

    # stator and magnet parameters of interior mounted magnet machine
    r = get_stator_dimensions(par)

    # magnet parameters
    r['magnet'] = get_interior_magnet_dimensions(
        r['i1'], r['w1'], r['kw'], r['psi1'], r['lfe'], r['Da2'],
        par)

    r['magnet'].pop('hM')
    _set_genpars(r, 2*par['p'])
    r['name'] = f"IPM-{r['poles']}"

    return r


def im(pnom, speed, p, **kwargs):
    """returns dimension of a IM machine"""
    par = dict(
        pnom=pnom, speed=speed, p=p)
    par.update(kwargs)
    _set_im_defaults(par)

    # stator and rotor parameters of induction machine
    try:
        slots = _im_stator_slots(p)
    except KeyError:
        slots = []
    r = get_stator_dimensions(par, slots=slots)
    # rotor parameters
    r['rotor'] = get_im_rotor_dimensions(
        par['cos_phi']*r['A'], r['Da2'], r['psi1'], r['lfe'], par)
    _set_genpars(r, 2*par['p'])
    r['name'] = f"IM-{r['poles']}"
    return r


def eesm(pnom, speed, p, **kwargs):
    """returns dimension of a EESM machine"""
    par = dict(
        pnom=pnom, speed=speed, p=p)
    par.update(kwargs)
    _set_sm_defaults(par)

    # stator and rotor parameters of ee synchronous machine
    r = get_stator_dimensions(par)

    # rotor parameters
    r['rotor'] = get_sm_rotor_dimensions(r['A'], r['psi1'],
                                         r['lfe'], r['Da2'], par)

    _set_genpars(r, 2*par['p'])
    r['name'] = f"EESM-{r['poles']}"
    return r


if __name__ == "__main__":

    pnom = 10e3
    speed = 4400/60
    p = 4
    udc = 600

    r = spm(pnom, speed, p, udc=udc)

    # print results
    print("Number of slots          {:10d}".format(r['stator']['num_slots']))
    print("Iron length            [mm] {:10.2f}".format(r['lfe']*1e3))
    print("Stator outer diameter  [mm] {:10.2f}".format(r['outer_diam']*1e3))
    print("Stator bore diameter   [mm] {:10.2f}".format(r['bore_diam']*1e3))
    print("Slot height            [mm] {:10.2f}".format(
        r['stator']['statorRotor3']['slot_height']*1e3))
#    print("     width             [mm] {:10.2f}".format(
#        r['stator']['statorRotor3']['slot_']*1e3))
    print("Current loading      [kA/m] {:10.2f}".format(r['A']*1e-3))
    print("Thermal load       [A²/mm³] {:10.2f}".format(r['AJ']*1e-9))
    print("Esson           [kVAmin/m³] {:10.2f}".format(r['ess']))
#    print("Windings                    {:7}".format(r['windings']))
    print("Windingfactor                {:10.3f}".format(r['kw']))
    print("Slots/Pole Phase            {:9.1f}".format(r['q']))
    print("Current                 [A] {:10.2f}".format(r['i1']))
    print("Iron Losses             [W] {:10.2f}".format(r['pfe']))
    print("Cur Losses              [W] {:10.2f}".format(r['pcu']))
#    json.dump(r, sys.stdout)
