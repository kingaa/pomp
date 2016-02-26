library(pomp)

pompExample(ou2)

set.seed(64857673L)

pdf(file="ou2-mif-fp.pdf")

p.truth <- coef(ou2)
guess2 <- guess1 <- p.truth
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.25*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 4*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

mif1 <- mif(ou2,Nmif=100,start=guess1,
            ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=5,x2.0=5,
              alpha.2=0.1,alpha.3=0.1
              ),
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.type="geometric",
            cooling.fraction=0.95^50,
            max.fail=100,
            method="fp"
            )

mif2 <- mif(ou2,Nmif=100,start=guess2,
            ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=5,x2.0=5,
              alpha.2=0.1,alpha.3=0.1
              ),
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.type="geometric",
            cooling.fraction=0.95^50,
            max.fail=100,
            method="fp"
            )

plot(c(mif1,mif2))

dev.off()
