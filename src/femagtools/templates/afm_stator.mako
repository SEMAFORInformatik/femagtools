--[[ Supported AFM Model Types
  "S1R1"      -- 1 stator, 1 rotor
  "S1R2"      -- 1 stator, 2 rotor, 1 half simulated
  "S1R2_all"  -- 1 stator, 2 rotor, all simulated
  "S2R1"      -- 2 stator, 1 rotor, 1 half simulated
  "S2R1_all"  -- 2 stator, 1 rotor, all simulated
--]]
m.slot_height     =  ${model['slot_height']*1e3}    -- Total slot height [mm]
m.slot_h1         =  ${model['slot_h1']*1e3}        -- Slot opening height [mm]
m.slot_h2         =  ${model['slot_h2']*1e3}        -- Distance to radius [mm]
m.st_yoke_height  =  ${model.get('yoke_height', 0)*1e3} -- Yoke height [mm]
m.slot_width      =  ${model['slot_width']*1e3}     -- Slot width [mm]
m.slot_open_width =  ${model['slot_open_width']*1e3} -- Slot opening width [mm]
m.slot_r1         =  ${model['slot_r1']*1e3}         -- Slot radius [mm]
m.slot_r2         =  ${model['slot_r2']*1e3}         -- Slot radius [mm]
m.mcvkey_yoke     =  mcvkey_yoke
m.nodedist        =  ${model.get('m.nodedist', 1)}  -- Node distance

  -- overdeterminate slot_h2
  if (math.abs(m.slot_h2 - m.slot_r1 - m.slot_h1) > 0) then
     m.slot_h2 = m.slot_h1 + m.slot_r1
  end
  -- overdeterminate slot_open_width
  if (m.slot_r1 > 0 and math.abs(m.slot_width - 2*m.slot_r1 - m.slot_open_width) > 0) then
     m.slot_open_width = m.slot_width - 2*m.slot_r1
  end

  if (m.model_type == "S1R2") or  (m.model_type == "S1R2_all") then
    if (m.st_yoke_height > 0) then
      m.st_yoke_height = 0
      print("Warning: m.st_yoke_height for this type must be 0. Set to 0!")
    end
  end

  if m.st_yoke_height > 0 then
    sh = m.slot_height
  else
    sh = m.slot_height/2
  end

  fml = require("fml")

  P0 = fml.Point:Create(0,0)

  Lh1 = fml.Line:Create(P0,0)
  Lh2 = fml.Line:Parallel(Lh1,-m.slot_h1)
  Lh3 = fml.Line:Parallel(Lh1,-m.slot_h2)
  Lh4 = fml.Line:Parallel(Lh1,-(sh-m.slot_r2))
  Lh5 = fml.Line:Parallel(Lh1,-sh)
  Lh6 = fml.Line:Parallel(Lh1,-(sh+m.st_yoke_height))

  Lv1 = fml.Line:Create(P0,90)
  m.width_bz = m.pole_width*m.num_poles/m.tot_num_slot
  m.tooth_width = m.width_bz-m.slot_width
  Lv2 = fml.Line:Parallel(Lv1,m.tooth_width/2)
  Lv3 = fml.Line:Parallel(Lv1,m.tooth_width/2+m.slot_r1)
  Lv4 = fml.Line:Parallel(Lv1,m.tooth_width/2+m.slot_r2)
  Lv5 = fml.Line:Parallel(Lv1,m.tooth_width/2+m.slot_width/2-m.slot_open_width/2)
  Lv6 = fml.Line:Parallel(Lv1,m.tooth_width/2+m.slot_width/2)

  P11 = P0
  P16 = fml.Point:Intersection(Lh6,Lv1)
  P23 = fml.Point:Intersection(Lh3,Lv2)
  P24 = fml.Point:Intersection(Lh4,Lv2)
  P25 = fml.Point:Intersection(Lh5,Lv2)
  P26 = fml.Point:Intersection(Lh6,Lv2)
  P34 = fml.Point:Intersection(Lh4,Lv3)
  P33 = fml.Point:Intersection(Lh3,Lv3)
  P44 = fml.Point:Intersection(Lh4,Lv4)
  P45 = fml.Point:Intersection(Lh5,Lv4)
  P51 = fml.Point:Intersection(Lh1,Lv5)
  P52 = fml.Point:Intersection(Lh2,Lv5)
  P61 = fml.Point:Intersection(Lh1,Lv6)
  P62 = fml.Point:Intersection(Lh2,Lv6)
  P65 = fml.Point:Intersection(Lh5,Lv6)
  P66 = fml.Point:Intersection(Lh6,Lv6)

  if (m.slot_r1 > 0) then
    C1 = fml.Circle:Create(P33,m.slot_r1)
    P32 = fml.Point:Tangent(P52,C1,1)
  end

  model_size = m.pole_width*m.num_poles*m.num_slots/m.tot_num_slot
  blow_up_wind(0,-model_size/2, model_size,model_size/2)
  ndt(m.nodedist)

  nc_line(P11.x,P11.y, P51.x,P51.y, 0)
  nc_line_cont(P52.x, P52.y, 0)
  if P32 ~= nil and m.slot_r1 > 0 then
    nc_line_cont(P32.x, P32.y, 0)
    nc_circle_m(P23.x,P23.y, P32.x,P32.y, P33.x,P33.y, m.slot_r1, 0)
  else
    nc_line_cont(P33.x, P33.y, 0)
  end
  if m.slot_r2 > 0 then
    nc_circle_m(P45.x,P45.y, P24.x,P24.y, P44.x,P44.y, m.slot_r2, 0)
    nc_line(P23.x, P23.y, P24.x, P24.y, 0)
  else
    if P32 ~= nil then
      nc_line(P23.x, P23.y, P45.x, P45.y, 0)
    else
      nc_line(P33.x, P33.y, P45.x, P45.y, 0)
    end
  end

  nc_line(P45.x, P45.y, P65.x, P65.y, 0)
  nc_line_cont(P66.x, P66.y, 0)
  nc_line_cont(P16.x, P16.y, 0)
  nc_line_cont(P11.x, P11.y, 0)

  x, y = (P11.x+P26.x)/2, (P11.y + P16.y)/2
  create_mesh_se(x, y)
  def_new_sreg(x,y, "stfe", "skyblue")

  nc_line(P51.x, P51.y, P61.x, P61.y, 0)
  nc_line_cont(P62.x, P62.y, 0)
  nc_line_cont(P52.x, P52.y, 0)
  if m.slot_h1 > 0 then
    create_mesh_se((P51.x+P61.x)/2, (P61.y + P62.y)/2)
  end
  nc_line(P62.x, P62.y, P65.x, P65.y, 0)
  m.xcoil_1, m.ycoil_1 = (P26.x+P62.x)/2, (P62.y + P65.y)/2
  create_mesh_se(m.xcoil_1, m.ycoil_1)

  x1,y1 = P61.x,P61.y
  x2,y2 = P66.x,P66.y
  mirror_nodechains(x1,y1, x2,y2)

  x1,y1 = P16.x,P16.y
  x2,y2 = P11.x,P11.y
  x3,y3 = 2*P61.x,P61.y
  x4,y4 = 2*P61.x,P66.y
  translate_copy_nodechains(x1,y1, x2,y2, x3,y3, x4,y4,
                  m.num_sl_gen-1)

  x,y = (P11.x+P26.x)/2,(P11.y+P26.y)/2
  if (m.mcvkey_yoke == 'dummy') then
    def_mat_fm(x,y, 'blue', 1000, 100)
  else
    def_mat_fm_nlin(x,y, 'blue', m.mcvkey_yoke, 100, 0)
  end

  m.middle_line  = ${model.get('middle_line', 0)}  -- vertical slot parts if 1
     for n = 1, m.num_sl_gen do
       s = 'S'..n
       x, y = (2*n-1)*m.width_bz/2 -(m.width_bz/2-m.xcoil_1), m.ycoil_1
       def_new_sreg(x, y, s, 'white')
       x = x + 2*(m.width_bz/2-m.xcoil_1)
       if m.middle_line == 0 then
         add_to_sreg(x, y, s)
       else
         def_new_sreg(x, y, s..'2', 'white')
       end
     end

    if m.middle_line == 1 then
      m.xcoil_2 = m.width_bz/2 + (m.width_bz/2-m.xcoil_1)
      m.ycoil_2 = m.ycoil_1
    end
