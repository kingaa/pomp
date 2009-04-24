require(pomp)

euler.sir <- simulate(
                      pomp(
                           times=seq(1/52,4,by=1/52),
                           data=rbind(measles=numeric(52*4)),
                           t0=0,
                           tcovar=seq(0,25,by=1/52),
                           covar=matrix(
                             periodic.bspline.basis(seq(0,25,by=1/52),nbasis=3,period=1,degree=3),
                             ncol=3,
                             dimnames=list(NULL,paste("seas",1:3,sep=''))
                             ),
                           delta.t=1/52/20,
                           statenames=c("S","I","R","cases","W","B","dW"),
                           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop","rho"),
                           covarnames=c("seas1"),
                           zeronames=c("cases"),
                           comp.names=c("S","I","R"),
                           step.fun="sir_euler_simulator",
                           rprocess=euler.simulate,
                           dens.fun="sir_euler_density",
                           dprocess=euler.density,
                           skeleton.vectorfield="sir_ODE",
                           rmeasure="binom_rmeasure",
                           dmeasure="binom_dmeasure",
                           PACKAGE="pomp",
                           initializer=function(params, t0, comp.names, ...){
                             p <- exp(params)
                             snames <- c(
                                         "S","I","R","cases","W","B",
                                         "SI","SD","IR","ID","RD","dW"
                                         )
                             fracs <- p[paste(comp.names,"0",sep=".")]
                             x0 <- numeric(length(snames))
                             names(x0) <- snames
                             x0[comp.names] <- round(p['pop']*fracs/sum(fracs))
                             x0
                           }
                           ),
                      params=log(
                        c(
                          gamma=26,mu=0.02,iota=0.01,
                          beta1=1200,beta2=1800,beta3=600,
                          beta.sd=1e-3,
                          pop=2.1e6,
                          rho=0.6,
                          S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
                          )
                        ),
                      nsim=1,
                      seed=329348545L
                      )
euler.sir <- euler.sir[[1]]
