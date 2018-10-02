options(digits=3)
png(filename="kalman2-%02d.png",res=100)

library(pomp2)

set.seed(1968726372)

gompertz() -> po

po %>%
  window(start=1,end=30) %>%
  as.data.frame() %>%
  subset(select=-X) -> dat

try(dat %>% enkf())

dat %>%
  enkf(
    times="time",t0=0,Np=100,
    params=c(r=0.1,K=150),
    rinit=function(K, ...) {
      c(x=K)
    },
    rprocess=discrete_time(
      function (x, r, K, ...) {
        e <- rnorm(n=1,mean=0,sd=0.1)
        c(x=x*exp(r*(1-x/K))+e)
      }
    ),
    R=2,
    h=function(x) 0.01*x
  ) %>% plot()

try(dat %>% eakf())

dat %>%
  eakf(
    times="time",t0=0,Np=100,
    params=c(r=0.1,K=100),
    rinit=function(K, ...) {
      c(x=rlnorm(n=1,meanlog=log(K),sdlog=0.01))
    },
    rprocess=discrete_time(
      function (x, r, K, ...) {
        c(x=x*exp(r*(1-x/K)))
      }
    ),
    R=0.01,
    C=0.01
  ) %>% plot()

try(po %>% enkf(rprocess=NULL))
try(po %>% eakf(rprocess=NULL))

dev.off()
