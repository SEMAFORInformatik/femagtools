"""Creating plots

"""
from .machine import machine
from .fluxdens import airgap, airgap_fft
from .bch import torque, torque_fft, force, force_fft, \
    fluxdens_surface, winding_current, winding_flux, \
    voltage, voltage_fft, \
    pmrelsim, multcal, fasttorque, cogging, transientsc, \
    transientsc_demag, i1beta_torque, i1beta_ld, i1beta_lq, i1beta_psid, \
    i1beta_psiq, i1beta_psim, i1beta_up, \
    idq_ld, idq_lq, idq_psid, idq_psim, idq_psiq, idq_torque
from .char import mtpa, mtpv, characteristics, efficiency_map, losses_map
from .forcedens import forcedens, forcedens_surface, forcedens_contour, forcedens_fft
from .nc import spel, mesh, demag, demag_pos, \
    flux_density, flux_density_eccentricity, \
    airgap_flux_density_pos, loss_density
from .mcv import mcv_hbj, mcv_muer, felosses
from .phasor import i1beta_phasor, iqd_phasor, phasor
from .wdg import mmf, mmf_fft, zoneplan, winding_factors, winding, currdist
from .fieldlines import fieldlines
