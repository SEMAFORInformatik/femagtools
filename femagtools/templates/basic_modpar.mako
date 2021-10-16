% if hasattr(model, 'num_agnodes'):
num_agnodes = ${model.num_agnodes}
agndst = 360/num_agnodes
% endif
% if hasattr(model, 'stator'):
% if hasattr(model, 'airgap'):
% if isinstance(model.get(['airgap']), list):
<%
ag = '{' + ','.join([str(x*1e3) for x in model.get(['airgap'])]) +'}'
%>
ag = ${ag}
% else:
ag  = ${model.get(['airgap'])*1e3}
% endif
% else:
ag = 0
% endif
% if model.move_action == 0:
% if model.external_rotor:
dy2 = ${model.get(['outer_diam'])*1e3}
% if isinstance(model.get(['bore_diam']), list):
<%
da2 = '{' + ','.join([str(x*1e3) for x in model.get(['bore_diam'])]) +'}'
%>
da2 = ${da2}
da1 = {}
for i=1, table.getn(ag) do
  da1[i] = da2[i] - 2*ag[i]
end
% else:  # bore_diam is scalar
da2 = ${model.get(['bore_diam'])*1e3}
da1 = da2 - 2*ag 
% endif
dy1 = ${model.get(['inner_diam'])*1e3}
% else: # internal rotor
dy1 = ${model.get(['outer_diam'])*1e3}
% if isinstance(model.get(['bore_diam']), list):
<%
da1 = '{' + ','.join([str(x*1e3) for x in model.get(['bore_diam'])]) +'}'
%>
da1 = ${da1}
da2 = {}
for i=1, table.getn(ag) do
  da2[i] = da1[i] - 2*ag[i]
end
% else: # bore_diam is scalar
da1 = ${model.get(['bore_diam'])*1e3}
da2 = da1 - 2*ag 
% endif
% if hasattr(model, 'shaft_diam'):
dy2 = ${model.get(['shaft_diam'])*1e3}
dsh = ${model.get(['inner_diam'])*1e3}
% else:
dy2 = ${model.get(['inner_diam'])*1e3}
% endif
% endif
% endif
% endif

% if hasattr(model, 'stator'):
% if 'num_slots' in model.stator: # DC motor with CMM
m.tot_num_slot    =   ${int(model.get(['stator','num_slots']))}
m.num_sl_gen      =   ${int(model.get(['stator','num_slots_gen']))}
% elif hasattr(model, 'rotor'):
m.tot_num_slot    =   ${int(model.get(['rotor','num_slots']))}
m.num_sl_gen      =   ${int(model.get(['rotor','num_slots_gen']))}
% endif
m.num_poles       =   ${int(model.get(['poles']))}
m.num_pol_pair    =   m.num_poles/2
m.num_slots       =   m.num_sl_gen 
m.npols_gen       =   m.num_poles * m.num_sl_gen / m.tot_num_slot
m.tot_num_sl      =   m.tot_num_slot
% if model.move_action == 0:
% if isinstance(model.get(['bore_diam']), list):
m.fc_radius       =   (da1[2]/2-ag[2]/2) -- Radius airgap (extern)
m.fc_radius1      =   (da1[1]/2-ag[1]/2) -- Radius airgap (intern?)
m.fc_radius2 = 	m.fc_radius1
% else:
m.fc_radius       =   (da1+da2)/4
m.fc_radius1      =   m.fc_radius
m.sl_radius       =   m.fc_radius      -- radius of sliding area 
% endif
% endif
m.arm_length      =   ${model.get(['lfe'])*1e3}
pre_models("basic_modpar")
% endif
