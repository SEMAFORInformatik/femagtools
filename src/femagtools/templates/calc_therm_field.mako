-- Stator aussen
-- Waermeuebergange einer längsangestroemten Platte
-- 10m/s => 57.0W/m2.K
-- 20m/s => 95.3W/m2.K
-- 30m/s => 130.0W/m2.K
dy1 = 269.24
da1 = 161.92
dy2 = 110.64
heat_trans = 57.0
area_factor = 10
rau = (dy1+0.5)/2.
xau,yau = pr2c(rau,math.pi/num_slots) -- Luft aussen gegen Umgebung
def_heat_transfer(xau,yau,skyblue,heat_trans,area_factor)

rau = (dy2-0.5)/2.
area_factor = 1
xau,yau = pr2c(rau,math.pi/num_slots)        -- Luft innen gegen Umgebung
def_heat_transfer(xau,yau,skyblue,heat_trans/5,area_factor)

---------------------------------------------
-- set boundary conditions -------------------
del_bcond()

beta = 2*math.pi*slots_gen/num_slots
h = 3.8
r1 = dy1/2 + h
r2 = dy2/2 - h
x0, y0 = pr2c(r1, beta)
x1, y1 = pr2c(r2, beta)
x2, y2 = r2, 0
x3, y3 = r1, 0
def_bcond_tp(x0, y0, x1, y1, x2, y2, x3, y3, 4)
def_bcond_to(x1, y1, x2, y2)
def_bcond_to(x3, y3, x0, y0)
--[[
ag  = 0.75
da2 = da1 - 2*ag
r1 = da2/2+ag/2    -- Inner airgap nodechain
r2 = da2/2+ag      -- Outer airgap nodechain

x0, y0 = r1, 0 -- pr2c(r1, beta)
x1, y1 = pr2c(r1, beta)

def_nccond_to(x0,y0, x1,y1)

x0, y0 = r2, 0
x1, y1 = pr2c(r2, beta)
def_nccond_to(r2,0, 0,r2)
--]]
---------------------------------------------
-- Verluste import losses
---------------------------------------------
import_losses_from_femag_dc()
color_gradation_th(0,0,tot,Losses,0,0,model.."_losses.svg")

---------------------------------------------
-- Temperatur berechnen ---------------------
---------------------------------------------
calc_therm_field()
color_gradation_th(0,0,tot,Temp,0.0,0,model.."_temp.svg")
