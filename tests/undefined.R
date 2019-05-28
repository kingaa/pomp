options(digits=3)

library(pomp)

pomp:::undefined()
pomp:::undefined("harry")
pomp:::undefined(function(x)x)
pomp:::undefined(pomp:::pomp_fun(function(x)x))
pomp:::undefined(pomp:::pomp_fun())
pomp:::undefined(parameter_trans())
pomp:::undefined(pomp:::rproc_plugin())
pomp:::undefined(pomp:::skel_plugin())
