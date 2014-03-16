if (Sys.getenv("POMP_FULL_TESTS")=="yes") {

  library(pomp)

  pompExample(ou2)

  pdf(file="ou2-mif2.pdf")

  p.truth <- coef(ou2)
  guess2 <- guess1 <- p.truth
  guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.9*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
  guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.2*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

  set.seed(64857673L)
  mif1a <- mif(ou2,Nmif=100,start=guess1,
               pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=.5,x2.0=.5,
                 alpha.2=0.1,alpha.3=0.1),
               transform=F,
               Np=1000,
               var.factor=1,
               ic.lag=10,
               cooling.type="hyperbolic",
               cooling.fraction=0.05,
               method="mif2",
               tol=1e-8
               )

  mif2a <- mif(ou2,Nmif=100,start=guess1,
               pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=0.5,x2.0=.5,
                 alpha.2=0.1,alpha.3=0.1),
               transform=F,
               Np=1000,
               var.factor=1,
               ic.lag=10,
               cooling.type="geometric",
               cooling.fraction=0.95^50,
               max.fail=100,
               method="mif",
               tol=1e-8
               )  

  compare.mif(list(mif1a,mif2a))

  set.seed(64857673L)
  mif1b <- mif(ou2,Nmif=50,start=guess1,
               pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=.5,x2.0=.5,
                 alpha.2=0.1,alpha.3=0.1),
               transform=F,
               Np=1000,
               var.factor=1,
               ic.lag=10,
               cooling.type="hyperbolic",
               cooling.fraction=0.05,
               method="mif2"
               )
  mif1b <- continue(mif1b,Nmif=50)

  mif2b <- mif(ou2,Nmif=50,start=guess1,
               pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=0.5,x2.0=.5,
                 alpha.2=0.1,alpha.3=0.1),
               transform=F,
               Np=1000,
               var.factor=1,
               ic.lag=10,
               cooling.whatsit=200,
               cooling.type="geometric",
               cooling.factor=0.95,
               max.fail=100,
               method="mif"
               )  
  mif2b <- continue(mif2b,Nmif=50)

  mif2c <- mif(ou2,Nmif=50,start=guess1,
               pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=0.5,x2.0=.5,
                 alpha.2=0.1,alpha.3=0.1),
               transform=F,
               Np=1000,
               var.factor=1,
               cooling.type="hyperbolic",
               cooling.fraction=0.05,
               max.fail=100,
               method="mif2"
               )  
  mif2c <- continue(mif2c,Nmif=50)

  compare.mif(list(mif1b,mif2b))

  compare.mif(list(mif1a,mif1b))
  compare.mif(list(mif2a,mif2b))

  compare.mif(list(mif1b,mif2c))

  dev.off()

}
