

ag  = ${model.get(['airgap'])*1e3}
% if model.external_rotor:
dy2 = ${model.get(['outer_diam'])*1e3}
da2 = ${model.get(['bore_diam'])*1e3}
dy1 = ${model.get(['inner_diam'])*1e3}
da1 = da2 - 2*ag 
% else:
dy1 = ${model.get(['outer_diam'])*1e3}
da1 = ${model.get(['bore_diam'])*1e3}
dy2 = ${model.get(['inner_diam'])*1e3}
da2 = da1 - 2*ag 
% endif

% if hasattr(model, 'stator'):
m.tot_num_slot    =   ${int(model.get(['stator','num_slots']))}
m.num_sl_gen      =   ${int(model.get(['stator','num_slots_gen']))}
m.num_poles       =   ${int(model.get(['poles']))}
m.num_pol_pair    =    m.num_poles/2
m.num_slots       =   m.num_sl_gen 
m.npols_gen       =   m.num_poles * m.num_sl_gen / m.tot_num_slot
m.tot_num_sl      =   m.tot_num_slot
% if model.move_action == 0:
m.fc_radius       =   (da1+da2)/4
m.fc_radius1      =   m.fc_radius
% endif
m.arm_length      =   ${model.get(['lfe'])*1e3}
pre_models("basic_modpar")
% endif

m.airgap          =     2*ag/3
