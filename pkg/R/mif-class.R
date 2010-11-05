## define the mif class
setClass(
         'mif',
         contains='pfilterd.pomp',
         representation=representation(
           ivps = 'character',
           pars = 'character',
           Nmif = 'integer',
           particles = 'function',
           var.factor='numeric',
           ic.lag='integer',
           cooling.factor='numeric',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )
