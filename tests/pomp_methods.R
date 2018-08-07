options(digits=3)

library(pomp)
library(magrittr)

set.seed(571163577)

pompExample(ou2,envir=NULL) -> ou2
ou2[[1]] -> po
po1 <- simulate(po)

as(po,"data.frame") %>% head()
as.data.frame(po1) %>% head()

obs(po)[,1:3]
obs(po,"y2")[,1:3]
try(obs(po,c("y2","z")))

states(po)
states(po1,"x1")[,1:3]
try(states(po1,"z"))
states(po1)[,1:3]

time(po)[1:3]
time(po,t0=TRUE)[1:3]

time(po) <- 1:10
try(time(po) <- c("bob","nancy"))
time(po1,t0=TRUE) <- 0:10
try(time(po) <- 10:0)
try(time(po,t0=TRUE) <- c(4,1:10))

window(po,end=5)
window(po,start=5)
window(po,start=5,end=10)
try(window(po,start=5,end=3))
try(window(po,start=NA,end=3))
try(window(po,start=1,end=NULL))

timezero(po)
timezero(po) <- -3
try(timezero(po) <- NA)
try(timezero(po) <- c(1,2,3))
try(timezero(po) <- 20)

coef(po)
coef(po,c("alpha.3","tau"))
try(coef(po,c("alpha.3","z")))

coef(po,"alpha.3") <- 4
coef(po,"z") <- 9
coef(po)
coef(po) <- NULL
coef(po)
coef(po) <- list(a=3,b=12)

pompExample(gompertz)
gompertz -> po

coef(po)
coef(po,transform=TRUE,pars=c("r","K"))
coef(po,"sigma",transform=TRUE) <- 0
coef(po)
coef(po,c("r","K")) <- c(a=1,b=2)
coef(po,transform=TRUE) <- c(r=1,K=1)
coef(po) <- NULL
try(coef(po,transform=FALSE) <- c(5,3))
try(coef(po,transform=TRUE) <- c(5,3))
coef(po,transform=TRUE) <- c(r=1,K=1)
coef(po)
po %>%
  window(start=5,end=20) %>%
  pomp(covar=data.frame(time=0:20,q=0:20),tcovar="time",
  larry=3L) -> po1
as(po1,"data.frame") %>% head()
