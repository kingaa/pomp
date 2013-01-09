library(pomp)

data(ou2)

set.seed(1285370209L)

pdf(file="aou2-mif2-mif.pdf")

p.truth <- coef(ou2)
guess2 <- guess1 <- p.truth
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.9*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.2*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

set.seed(64857673L)
mif1 <- mif(ou2,Nmif=100,start=guess1,
            pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=.5,x2.0=.5,
              alpha.2=0.1,alpha.3=0.1),
            transform=F,
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.factor=0.95,
            cooling.fraction=0.05,
            option="mif2",
            .ndone=0,tol=10^-8,
            verbose=F
            )
##set.seed(64857673L)
mif2 <- mif(ou2,Nmif=100,start=guess1,
            pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=0.5,x2.0=.5,
              alpha.2=0.1,alpha.3=0.1),
            transform=F,
            cooling.scalar=430,
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.factor=0.95,
            cooling.fraction=0.5,
            max.fail=100,
            option="mif",
            .ndone=0,tol=10^-8,
            verbose=F
            )  
compare.mif(list(mif1,mif2))

dev.off()
