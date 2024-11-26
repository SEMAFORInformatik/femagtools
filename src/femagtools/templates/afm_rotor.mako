-- Create model of the rotor of the axial flow machine

m.yoke_height  = ${model.get('yoke_height', 0)*1e3}  -- Yoke height [mm]
m.magn_height     = ${model['magn_height']*1e3}  -- Magnet height [mm]
m.magn_width      = ${model['rel_magn_width']*1e2}  -- Magnet width: >0:[%]; <0: [mm]
m.spoke_width     = ${model.get('spoke_width', 0)*1e3}  -- Spoke width: >0 magn_width will be calculated
m.gap_ma_yoke     = 0
m.nodedist        = ${model.get('nodedist', 1)}

  if (m.model_type == "S2R1") or  (m.model_type == "S2R1_all") then
    if (m.yoke_height > 0) then
      m.yoke_height = 0
      print("Warning: m.yoke_height for this typ must be 0. Set to 0!")
    end
  end

  if (m.yoke_height > 0) then
    mh = m.magn_height
  else
    mh = m.magn_height/2
  end

  if m.spoke_width > 0 then
    m.magn_width = -(m.pole_width-m.spoke_width)
  end
  if m.magn_width > 0 then
    mw = m.magn_width/100*m.pole_width
  else
    mw = -m.magn_width
  end

  fml = require("fml")

  P0 = fml.Point:Create(0, -ag)

  Lh1 = fml.Line:Create(P0,0)
  Lh2 = fml.Line:Parallel(Lh1,mh)
  Lh3 = fml.Line:Parallel(Lh2,m.yoke_height)

  Lv1 = fml.Line:Create(P0,90)
  Lv2 = fml.Line:Parallel(Lv1,(m.pole_width-mw)/2)
  Lv3 = fml.Line:Parallel(Lv1,(m.pole_width+mw)/2)
  Lv4 = fml.Line:Parallel(Lv1,m.pole_width)

  P11 = P0
  P12 = fml.Point:Intersection(Lh2,Lv1)
  P13 = fml.Point:Intersection(Lh3,Lv1)
  P21 = fml.Point:Intersection(Lh1,Lv2)
  P22 = fml.Point:Intersection(Lh2,Lv2)
  P31 = fml.Point:Intersection(Lh1,Lv3)
  P32 = fml.Point:Intersection(Lh2,Lv3)
  P41 = fml.Point:Intersection(Lh1,Lv4)
  P42 = fml.Point:Intersection(Lh2,Lv4)
  P43 = fml.Point:Intersection(Lh3,Lv4)

  ndt(m.nodedist)

  nc_line(P11.x,P11.y, P21.x,P21.y, 0)
  nc_line_cont(P31.x,P31.y, 0)
  nc_line_cont(P41.x,P41.y, 0)

  nc_line(P12.x,P12.y, P22.x,P22.y, 0)
  nc_line_cont(P32.x,P32.y, 0)
  nc_line_cont(P42.x,P42.y, 0)

  nc_line(P11.x,P11.y, P12.x,P12.y, 0)
  nc_line(P21.x,P21.y, P22.x,P22.y, 0)
  nc_line(P31.x,P31.y, P32.x,P32.y, 0)
  nc_line(P41.x,P41.y, P42.x,P42.y, 0)

  create_mesh_se((P12.x+P22.x)/2, (P11.y + P12.y)/2)
  create_mesh_se((P22.x+P32.x)/2, (P21.y + P22.y)/2)
  create_mesh_se((P32.x+P42.x)/2, (P31.y + P32.y)/2)

  if (m.yoke_height > 0) then
    nc_line(P13.x,P13.y, P43.x,P43.y, 0)
    nc_line(P12.x,P12.y, P13.x,P13.y, 0)
    nc_line(P42.x,P42.y, P43.x,P43.y, 0)
    x,y = (P12.x+P42.x)/2,(P42.y+P43.y)/2
    create_mesh_se(x,y)
    def_new_sreg(x,y, "rofe", "skyblue")
    if (mcvkey_yoke == 'dummy') then
      def_mat_fm(x,y, 'blue', 1000, 100)
    else
      def_mat_fm_nlin(x,y, 'blue', mcvkey_yoke, 100, 0)
    end
  end

  translate_copy_nodechains(P11.x,P11.y, P13.x,P13.y, P43.x,P43.y, P41.x,P41.y, m.npols_gen-1)

  x,y = (P22.x+P32.x)/2,(P21.y+P22.y)/2
  if m.remanenc == nil then
    m.remanenc = 1.2
  end
  if m.relperm == nil then
    m.relperm = 1.05
  end
  if m.rlen == nil then
    m.rlen = 100
  end
  for i = 1,m.npols_gen do
    if i%2 == 0 then
      def_mat_pm(x,y, 'red', m.remanenc, m.relperm, 90, "parallel", 0, m.rlen)
    else
      def_mat_pm(x,y, 'green', m.remanenc, m.relperm, 270, "parallel", 0, m.rlen)
    end
    x = x + m.pole_width
  end
