-- EXPERIMENTAL
-- Caution: this will likely be modified
--
local open = io.open

local function read_load_file(path)
    local file = open(path, "r") -- r read mode and b binary mode
    --if not file then return nil end
    local load = {}

    for line in io.lines(path) do
    local t, n, T, i1, beta = line:match("([^,]*),([^,]*),([^,]*),([^,]*),([^,]*)")
    load[#load+1] = { t=tonumber(t), n = tonumber(n), T = tonumber(T), i1 = tonumber(i1), beta = tonumber(beta) }
    end

    --file:close()
    return load
end

state_of_problem("therm_static")

beta = 360*m.npols_gen/m.num_poles
dy1 = ${model['outer_diam']*1e3}
dy2 = ${model['inner_diam']*1e3}
m.zeroangl = 0
-- heat transfer outside
heat_transfer_coefficient = ${model['heat_transfer_coefficient'][0]}
area_factor = 1
x,y = pd2c(dy1/2+.05, beta/2+m.zeroangl)
def_heat_transfer(x,y,'yellow',heat_transfer_coefficient, area_factor)
-- heat transfer inside
x,y = pd2c(dy2/2-0.05, beta/2+m.zeroangl)
heat_transfer_coefficient = ${model['heat_transfer_coefficient'][1]}
def_heat_transfer(x,y,'yellow',heat_transfer_coefficient, area_factor)

--
-- Torq Calc
--
m.move_action     =    0.0 -- rotate
m.skew_angle      =    0
m.nu_skew_steps   =    0
m.magn_temp       =   ${model['magn_temp']}
m.winding_temp   =    ${model['wind_temp']}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.fc_force_points   =  0.0 --    number move points in air gap
m.nu_force_pat   =        0  --   Number of extra force pathes < = 3
m.phi_start       =    0.0 --
m.range_phi       =    720./m.num_poles
m.nu_move_steps   =    30
m.num_par_wdgs    =    ${model['num_par_wdgs']}
m.period_frac     = 6

m.pocfilename    = model..'_'..m.num_poles..'p.poc'
xcoil, ycoil = ${model['temp_coords'][0]},${model['temp_coords'][1]}

local load = read_load_file("load.csv")

f = io.open('temperature.dat', 'w')
f:write(string.format("0.0 0.0\n"))

i = 1
for _, load_data in ipairs(load) do  -- use pairs or ipairs to iterate over tables
  if i == 1 then
    t1 = load_data.t
    i1 = load_data.i1
    beta = load_data.beta
    speed = load_data.n
    printf("\nTotal steps: %d\n", #load)
  else
    t2 = load_data.t
    end_time = t2 - t1
    i1 = load_data.i1
    beta = load_data.beta
    speed = load_data.n

    m.current   = i1 * math.sqrt(2)/m.num_par_wdgs
    m.angl_i_up = beta
    m.speed     = speed
    printf("Step %3d: %6.1fs %5.1f A %5.1f° %6.1f rpm: ",
            i-1, t1, i1, m.angl_i_up, m.speed)

    state_of_problem("mag_static")
    run_models("torq_calc")

    state_of_problem("therm_static")
    import_losses_from_femag_dc()

    if i == 2 then
      start_mode = temp_reset
    else
      start_mode = temp_load
    end
    time_step = end_time/5
    calc_therm_field_tstep(start_mode,time_step,end_time)

    export_calc_results("temp-"..string.format("%03d",i)..".vtu")
    temp = temperature_xy(xcoil, ycoil)
    printf("%4d, time: %6.1fs, current: %5.1fA, speed: %6.1frpm, Temp incr.: %2.2fK, time-step: %3.1fs\n", i-1, t2, i1, m.speed, temp, time_step)
    f:write(string.format("%g %g\n", t2, temp))
    t1 = t2
  end

  i = i+1

end
f:close()
save_model('close')
