## define the mif class
setClass(
         'mif',
         representation(
                        ivps = 'character',
                        pars = 'character',
                        Nmif = 'integer',
                        particles = 'function',
                        alg.pars = 'list',
                        random.walk.sd = 'numeric',
                        pred.mean = 'matrix',
                        pred.var = 'matrix',
                        filter.mean = 'matrix',
                        conv.rec = 'matrix',
                        eff.sample.size = 'numeric',
                        cond.loglik = 'numeric',
                        loglik = 'numeric'
                        ),
         contains='pomp'
         )
