png(filename="baddm-%02d.png",res=100)
options(digits=2)

library(pomp)

set.seed(398585L)

try(
  gompertz() %>%
    pfilter(Np=10,dmeasure=function(...,log) Inf)
)

try(
  gompertz() %>%
    mif2(
      Np=10,Nmif=2,
      dmeasure=function(...,log) NaN,
      rw.sd=rw.sd(r=0.01),
      cooling.fraction.50=0.1
    )
)

try(
  gompertz() %>%
    bsmc2(
      Np=10,
      rprior=function(...)c(r=runif(n=1,min=0,max=1)),
      dmeasure=function(...,log) NA_real_,
      partrans=NULL
    )
)

dev.off()
