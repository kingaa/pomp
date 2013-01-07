## define the mif class
setClass(
         'mif',
         contains='pfilterd.pomp',
         representation=representation(
           transform = "logical",
           ivps = 'character',
           pars = 'character',
           Nmif = 'integer',
           option='character',
           cooling.fraction='numeric',
           particles = 'function',
           var.factor='numeric',
           ic.lag='integer',
           cooling.factor='numeric',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )
