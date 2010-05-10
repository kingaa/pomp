require(pomp)

simulate(
         pomp( 
              times=seq(1,100),
              data=rbind(
                y1=rep(0,100),
                y2=rep(0,100)
                ),
              t0=0,
              rprocess = function (xstart, times, params, paramnames, ...) {
                nvar <- nrow(xstart)
                npar <- nrow(params)
                nrep <- ncol(xstart)
                ntimes <- length(times)
                ## get indices of the various parameters in the 'params' matrix
                ## C uses zero-based indexing!
                parindex <- match(paramnames,rownames(params))-1
                array(
                      .C("ou2_adv",
                         X = double(nvar*nrep*ntimes),
                         xstart = as.double(xstart),
                         par = as.double(params),
                         times = as.double(times),
                         n = as.integer(c(nvar,npar,nrep,ntimes)),
                         parindex = as.integer(parindex),
                         DUP = FALSE,
                         NAOK = TRUE,
                         PACKAGE = "pomp"
                         )$X,
                      dim=c(nvar,nrep,ntimes),
                      dimnames=list(rownames(xstart),NULL,NULL)
                      )
              },
              dprocess = function (x, times, params, log, paramnames, ...) {
                nvar <- nrow(x)
                npar <- nrow(params)
                nrep <- ncol(x)
                ntimes <- length(times)
                parindex <- match(paramnames,rownames(params))-1
                array(
                      .C("ou2_pdf",
                         d = double(nrep*(ntimes-1)),
                         X = as.double(x),
                         par = as.double(params),
                         times = as.double(times),
                         n = as.integer(c(nvar,npar,nrep,ntimes)),
                         parindex = as.integer(parindex),
                         give_log=as.integer(log),
                         DUP = FALSE,
                         NAOK = TRUE,
                         PACKAGE = "pomp"
                         )$d,
                      dim=c(nrep,ntimes-1)
                      )
              },
              dmeasure = "normal_dmeasure",
              rmeasure = "normal_rmeasure",
              skeleton.map = function (x, t, params, ...) {
                with(
                     as.list(c(x,params)),
                     {
                       c(
                         x1=alpha.1*x1+alpha.3*x2,
                         x2=alpha.2*x1+alpha.4*x2
                         )
                     }
                     )
              },
              paramnames = c(
                "alpha.1","alpha.2","alpha.3","alpha.4",
                "sigma.1","sigma.2","sigma.3",
                "tau"
                ),
              statenames = c("x1","x2")
              ),
         params=c(
           alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
           sigma.1=3, sigma.2=-0.5, sigma.3=2,
           tau=1, 
           x1.0=-3, x2.0=4
           ),
         nsim=1,
         seed=377456545L
         ) -> ou2
