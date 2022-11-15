
-- conductor properties
m.conduct = ${model.get('elconduct', 10e6)} -- el. conductivity S/m
m.conrelperm = ${model.get('relperm', 100)} -- rel permeability
m.contemp = ${model.get('contemp', 20)} -- conductor temperature °C
m.contecoef = ${model.get('tempcoef', 0)} -- conduct. temperature coeff 1/K
m.spconweight = ${model.get('spmaweight', 7.6)} -- mass density g/cm³
m.relconlength = ${model.get('relconlength', 100)} -- rel cond length %

pre_models('conduct-data')