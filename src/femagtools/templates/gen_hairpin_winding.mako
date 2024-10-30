-- hairpin winding template

slot_div_angle = 360/m.tot_num_slot

if m.xcoil_1 ~= nil and m.ycoil_1 ~= nil then 
  rcoil, phi_coil = c2pd(m.xcoil_1, m.ycoil_1)
else
  rcoil, phi_coil = da1/2 + m.slot_height/2, slot_div_angle/2
end 

-- delete existing mesh in the slot
for i = 1, m.num_sl_gen do
    xc, yc = pd2c(rcoil, phi_coil+(i-1)*slot_div_angle)
    delete_mesh_se(xc, yc)
end
-- wire params
wire_height = ${model['wire']["wire_height"]*1e3}
wire_width = ${model['wire']["wire_width"]*1e3}
wire_corner_r = 0.0
wire_gap = ${model['wire']["wire_separation"]*1e3}
TH1 = ${model['wire']["wire_th1"]*1e3} -- distance to the bore radius
nwires = ${model['wire']["num_layers"]}

fml = require "fml"

function line_xy(pnt1, pnt2, nnodes)
   -- create line
  nc_line(pnt1.x, pnt1.y, pnt2.x, pnt2.y, nnodes)
end

function circle_xy(pnt1, pnt2, cx, cy, nnodes)
   -- create circle
  nc_circle_m(pnt1.x, pnt1.y, pnt2.x, pnt2.y, cx, cy, nnodes)
end

function draw_rectangle_wire(dx, dy, dphi, wire_height, wire_width, r0, agndst)

    local x = {}
    local y = {}
    local x1 = {}
    local y1 = {}
    local hl = {}
    local hi = {}
    local hj = {}
    local hk = {}
    local xcp, ycp
    local tmpx, tmpy, tmpr, tmpp

    h1 = fml.Point:Create(dx, 0)
    h2 = fml.Point:Create(dx, -wire_width/2+dy)
    he = fml.Point:Create(dx+wire_height, 0)

    hl[1] = fml.Line:Create(h1, 90)
    hl[2] = fml.Line:Parallel(hl[1], wire_height)
    hl[3] = fml.Line:Create(h2, 0)
    hl[4] = fml.Line:Parallel(hl[3], -wire_width)

    hi[1] = fml.Point:Intersection(hl[1], hl[3])
    hi[2] = fml.Point:Intersection(hl[1], hl[4])
    hi[3] = fml.Point:Intersection(hl[2], hl[4])
    hi[4] = fml.Point:Intersection(hl[2], hl[3])
    -- with corner
    if r0 > 0 then
       s1 = {1, 1, 2, 2}
       s2 = {3, 4, 4, 3}
       s3 = {1, 4, 3, 2}
       ctr = 1

       for i =1, 4 do
          hj[ctr], hj[ctr+1], hj[ctr+2] = fml.Point:Rounding(hl[s1[i]], hl[s2[i]], r0, s3[i])
          ctr = ctr + 3
       end

       for j = 1, #hj do
          tmpr, tmpp = c2pd(hj[j].x, hj[j].y)
          tmpx, tmpy = pd2c(tmpr, tmpp+dphi)
          hk[j] = fml.Point:Create(tmpx, tmpy)
       end

       dst = agndst/2
       ah1 = fml.Point:Create((hk[12].x+hk[8].x)/2 , (hk[12].y+hk[8].y)/2)
       ah2 = fml.Point:Create((hk[6].x+hk[2].x)/2 , (hk[6].y+hk[2].y)/2)

       line_xy(hk[12], ah1, 0)
       line_xy(ah1, hk[8], 0)
       line_xy(hk[9], hk[5], 0)
       line_xy( hk[6], ah2, 0)

       line_xy(ah2, hk[2], 0)
       line_xy(hk[3], hk[11], 0)

       circle_xy(hk[8], hk[9], hk[7].x, hk[7].y, 6)
       circle_xy(hk[5], hk[6], hk[4].x, hk[4].y, 6)
       circle_xy(hk[2], hk[3], hk[1].x, hk[1].y, 6)
       circle_xy(hk[11], hk[12], hk[10].x, hk[10].y, 6)
       xcp = (hk[12].x + hk[6].x)/2
       ycp = (hk[12].y + hk[6].y)/2

    else
      xcp = 0
      ycp = 0

      for i = 1, 4 do
          x[i], y[i] = hi[i].x, hi[i].y
          tmpr, tmpp = c2pd(x[i], y[i])
          tmpx, tmpy = pd2c(tmpr, tmpp+dphi)
          x1[i], y1[i] = tmpx, tmpy
          xcp = xcp + tmpx/4
          ycp = ycp + tmpy/4
      end

      nc_line(x1[1], y1[1], x1[2], y1[2], 0)
      nc_line_cont(x1[3], y1[3], 0)
      nc_line_cont(x1[4], y1[4], 0)
      nc_line_cont(x1[1], y1[1], 0)

    end

    return xcp, ycp, ah1, ah2
end

-- hairpin winding
wire_xy = {}
x0, y0 = da1/2+TH1, 0
auh1 = {}
auh2 = {}

for j = 1, m.num_sl_gen do
   wire_xy[j] = {}
   for k = 1, nwires do
    wire_xy[j][k] = 0.0
   end
end

ndt(agndst*1.5)
for j = 1, m.num_sl_gen do
  for i = 1, nwires do
    dx, dy = x0+(wire_gap+wire_height)*(i-1), y0
    xcp, ycp, ah1, ah2 = draw_rectangle_wire(dx, dy, 360/m.tot_num_slot/2*(2*j-1), wire_height, wire_width, wire_corner_r, agndst)
    wire_xy[j][nwires - (i-1)] =  fml.Point:Create(xcp, ycp) -- large radius to smaller radius
    create_mesh_se(xcp, ycp)
  end

end
create_mesh()

-- create winding

--widfile = io.open("wid.fsl", 'r')
widfile = io.open("wid.txt", 'r')
nrows = 1
winding = {}

for i in widfile:lines() do
  winding[nrows] = {}
  ncols = 1
  for j in string.gmatch(i, "[+-]?%d+") do
    winding[nrows][ncols] = tonumber(j)
    ncols = ncols + 1
  end
  nrows = nrows + 1
end

--[[
winding table
key slot layer dir
]]--
dir = 'wi'
cols = {"green", "yellow", "cyan"}
wk = 1
for i = 1, #winding do
    windingskey = winding[i][1]
    slot_nr = winding[i][2]
    layer = winding[i][3]
    direction = winding[i][4]
    if direction < 0 then 
      dir = 'wo'
    else
      dir = 'wi'
    end 
    if i == 1 then
          if slot_nr <= m.num_sl_gen then
              wkey = def_new_wdg(wire_xy[slot_nr][layer].x, wire_xy[slot_nr][layer].y, cols[winding[i][1]], "Phase"..winding[i][1], 1, 0, 0, dir)
          end

    else
          if winding[i][1] == winding[i-1][1] and winding[i][1] == wk then
              if slot_nr <= m.num_sl_gen then
                  add_to_wdg (wire_xy[slot_nr][layer].x, wire_xy[slot_nr][layer].y, "wsamekey", dir, "wser")
              end
            
          else
              if slot_nr <= m.num_sl_gen then
                  wkey = def_new_wdg(wire_xy[slot_nr][layer].x, wire_xy[slot_nr][layer].y, cols[winding[i][1]], "Phase"..winding[i][1], 1, 0, 0, dir)
                  wk = wk + 1
              end
          end
   end
end
m.num_par_wdgs = ${model.get('num_par_wdgs', 1)}