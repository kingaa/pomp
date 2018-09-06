options(digits=3)

library(pomp2)

pomp2:::undefined()
pomp2:::undefined("harry")
pomp2:::undefined(function(x)x)
pomp2:::undefined(pomp2:::pomp_fun(function(x)x))
pomp2:::undefined(pomp2:::pomp_fun())
pomp2:::undefined(parameter_trans())
pomp2:::undefined(pomp2:::rproc_plugin())
pomp2:::undefined(pomp2:::skel_plugin())
