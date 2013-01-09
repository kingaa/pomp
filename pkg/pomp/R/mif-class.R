## define the mif class
setClass(
         'mif',
         contains='pfilterd.pomp',
         representation=representation(
           transform = "logical",
           ivps = 'character',
           pars = 'character',
           Nmif = 'integer',
           particles = 'function',
           var.factor='numeric',
           ic.lag='integer',
           cooling.factor='numeric',
           cooling.fraction='numeric',
           method='character',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )
