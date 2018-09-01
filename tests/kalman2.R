options(digits=3)

library(pomp)
library(magrittr)

set.seed(1968726372)

pompExample(gompertz)

gompertz %>%
  window(start=1,end=10) %>%
  as.data.frame() -> dat

try(dat %>% enkf())

dat %>%
  enkf(
    times="time",t0=0,Np=100,
    params=c(r=0.1,K=150),
    rinit=function(K, ...) {
      c(x=rlnorm(n=1,meanlog=log(K),sdlog=1))
    },
    rprocess=discrete_time(
      function (x, r, K, ...) {
        c(x=x*exp(r*(1-x/K)))
      }
    ),
    R=2,
    h=function(x) 0.01*x
  ) %>% plot()

try(dat %>% eakf())

dat %>%
  subset(time<10) %>%
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
