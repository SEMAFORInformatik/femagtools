  exit_on_error=false
  exit_on_end=false
  verbosity=2

model = '${model.get(['name'])}'
load_model(model)

  global_unit(mm)
  pickdist(0.001)
  cosys(polar)

<%include file="fe-contr.mako" />

