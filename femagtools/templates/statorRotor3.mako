% if model.get('is_rotor', False):
m.yoke_diam   = dy2
m.inside_diam = da2
m.nodedist        =   ${model.get('nodedist',1)}
Q2 = ${int(model['num_slots'])}
Q1 = m.tot_num_sl
m.tot_num_sl  =   Q2
% if model.get('num_slots_gen', 0):
m.num_sl_gen  =   ${model['num_slots_gen']}
% else:
m.num_sl_gen  =   Q2 * m.npols_gen/m.num_poles
% endif
% else:
m.yoke_diam   = dy1
if (type(da1) == "table") then
  m.inside_diam = da1[2]
else
  m.inside_diam = da1
end
m.wdg_location = -1 -- for gen_windings
% endif
m.slot_height = ${model['slot_height']*1e3}
m.slot_h1     = ${model['slot_h1']*1e3}
m.slot_h2     = ${model['slot_h2']*1e3}
m.slot_width  = ${model['slot_width']*1e3}
m.slot_r1     = ${model['slot_r1']*1e3}
m.slot_r2     = ${model['slot_r2']*1e3}
m.wedge_width1= ${model['wedge_width1']*1e3}
m.wedge_width2= ${model['wedge_width2']*1e3}
m.middle_line = ${model.get('middle_line',0)}
m.tooth_width = ${model['tooth_width']*1e3}
m.slot_top_sh = ${model['slot_top_sh']}

m.zeroangl    = ${model.get('zeroangle',0)}
m.rlength     = ${model.get('rlength',1)*100}

m.mcvkey_yoke = mcvkey_yoke

pre_models("STATOR_3")

% if model.get('is_rotor', False):

if dsh ~= nil then -- create shaft subregion
  del_bcond()
  alph0=m.zeroangl/180*math.pi
  rlen=100
  x = {}
  y = {}
  r = dsh/2
  minr = 5e-2
  if(dsh/dy2 < minr) then
    r = dy2/8
  end
  taus = 2*math.pi/Q2
  if(m.num_sl_gen < m.tot_num_sl) then
    x[1],y[1] = pr2c(r,alph0)
    x[2],y[2] = pr2c(dy2/2,alph0)
    ndt(5*agndst)
    nc_line(x[1], y[1], x[2], y[2], 0)
    for i=1, m.num_sl_gen do
      x[3],y[3] = pr2c(r,taus*i+alph0)
      x[4],y[4] = pr2c(dy2/2,taus*i+alph0)
      nc_line(x[4], y[4], x[3], y[3], 0)
      nc_line_cont(x[1], y[1], 0)
      x[1], y[1] = x[3], y[3]
      x0, y0 = pr2c((r+dy2/2)/2, (2*i-1)*taus/2+alph0)
      create_mesh_se(x0, y0)
      if(i==1) then
        def_new_sreg(x0, y0,"Shft",'lightgrey')
      else
        add_to_sreg(x0, y0)
      end
    end

    if(dsh/dy2 < minr) then
      if(dsh >1) then
        x1, y1 = pr2c(dsh/2, alph0)
        x2, y2 = pr2c(dsh/2, m.num_sl_gen*taus + alph0)
        n = math.max(math.pi*dsh/agndst/5, 16)*m.num_sl_gen/Q2
        nc_circle_m(x1, y1, x2, y2, 0.0, 0.0, n)
        x3, y3 = pr2c(r, alph0)
        nc_line(x1, y1, x3, y3, 0)
        x3, y3 = pr2c(r, m.num_sl_gen*taus + alph0)
        nc_line(x2, y2, x3, y3, 0)
        x0, y0 = pr2c((r+dsh/2)/2, taus+alph0)
      else
        x0, y0 = pr2c(dy2/10, taus+alph0)
        x1, y1 = pr2c(r, alph0)
        x2, y2 = pr2c(r, m.num_sl_gen*taus + alph0)
        nc_line(0, 0, x1, y1, 0)
        nc_line(0, 0, x2, y2, 0)
      end
      create_mesh_se(x0, y0)
      add_to_sreg(x0, y0)
      if(dsh >1) then
        def_bcond_vpo(x1, y1, x2, y2)
      end
    end
  else -- full model
      if(dsh >1) then
        x1, y1 = pr2c(dsh/2, alph0)
        x2, y2 = pr2c(dsh/2, m.num_sl_gen*taus/2 + alph0)
        n = math.max(math.pi*dsh/agndst/5, 16)*m.num_sl_gen/Q2/2
        nc_circle_m(x1, y1, x2, y2, 0.0, 0.0, n)
        nc_circle_m(x2, y2, x1, y1, 0.0, 0.0, n)
        x0, y0 = pr2c((dy2+dsh)/4, taus+alph0)
        create_mesh_se(x0, y0)
        def_bcond_vpo(x1, y1, x2, y2)
        def_bcond_vpo(x2, y2, x1, y1)
      else
        x0, y0 = pr2c((dy2+dsh)/4, taus)
        create_mesh_se(x0,y0)
      end
      def_new_sreg(x0, y0,"Shft",'white')
  end
  x0, y0 = pr2c((dy2+dsh)/4, taus/2+alph0)
  if mcvkey_shaft ~= nil and mcvkey_shaft ~= 'dummy' then
    def_mat_fm_nlin(x0, y0, "lightgrey", mcvkey_shaft, rlen)
  else
    def_mat_fm(x0, y0, "lightgrey", 1000, rlen)
  end
end

r=da2/2-m.slot_height/2
phi=math.pi/Q2
x,y=pr2c(r,phi)
def_new_subreg( x,y, 'Bar', violet )

for i=2, m.num_sl_gen do
  phi=phi+2*math.pi/Q2
  x,y=pr2c(r,phi)
  add_to_subreg( x, y )
 end

% else:
if mcvkey_teeth ~= nil then
  if m.inside_diam > m.yoke_diam then
     r = (m.inside_diam - m.slot_height)/2
  else
     r = (m.inside_diam + m.slot_height)/2
  end
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end
%endif