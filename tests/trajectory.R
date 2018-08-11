options(digits=3)
png(filename="trajectory-%02d.png",res=100)

library(pomp)
library(magrittr)

pomp(
  data=NULL,
  times=seq(0,60,by=0.1),t0=0,
  skeleton=vectorfield(Csnippet("
    DV = c*(V-pow(V,3)/3 - R + i);
    DR = (V + a - b*R)/c;"
  )),
  rinit=Csnippet("
    V = 1; R = 0;"
  ),
  statenames=c("V","R"),
  paramnames=c("c","i","a","b"),
  params=c(a=0.7,b=0.8,c=2,i=0.8)
) -> fhn

x <- array(c(0,1,1,2,1,1,0,-1),
  dim=c(2,2,2),
  dimnames=list(c("V","R"),NULL,NULL))
params <- rbind(a=c(0.7,0.5),b=c(0.8,0.5),c=c(2,5),i=c(0.8,0))
skeleton(fhn,x,t=c(0,3),params=params)
y <- trajectory(fhn,params=params,hmax=0.1)
invisible(y[,,599:601])
matplot(time(fhn),t(y["V",,]),type='l',lty=1)
plot(y[1,,],y[2,,],type='n')
points(y[1,1,],y[2,1,],pch='.',cex=3,col='black')
points(y[1,2,],y[2,2,],pch='.',cex=3,col='red')

## nonautonomous case
pomp(
  data=fhn,
  times=seq(0,100,by=0.01),t0=0,
  tcovar=seq(0,101,by=0.1),
  covar=cbind(i=sin(2*pi*seq(0,101,by=0.1))),
  skeleton=vectorfield(
    function(x,t,params,covars,...) {
      c(
        V=unname(params['c']*(x['V']-(x['V']^3)/3-x['R']+covars['i'])),
        R=unname((x['V']+params['a']-params['b']*x['R'])/params['c'])
      )
    }
  )
) -> fhn1

invisible(skeleton(fhn1,x,t=c(0,3),params=params))
y <- trajectory(fhn1,params=params,hmax=0.01)
invisible(y[,,199:201])
matplot(time(fhn1),t(y["V",,]),type='l',lty=1)
plot(y[1,,],y[2,,],type='n')
points(y[1,1,],y[2,1,],pch='.',cex=3,col='black')
points(y[1,2,],y[2,2,],pch='.',cex=3,col='red')

invisible(trajectory(fhn,times=c(1,5)))
try(trajectory(fhn,times=NULL))
try(trajectory(fhn,times=c(1,1,1)))
try(trajectory(fhn,t0=10))
try(trajectory(fhn,params=c(3,2,1)))
try(trajectory(fhn,params=matrix(c(3,2,1,5),2,2)))
try(trajectory(fhn,params=NULL))
try(trajectory(fhn,params=list(a=3,b=2)))
try(trajectory(fhn,maxsteps=-1))
try(trajectory(fhn,maxsteps=1,verbose=TRUE) -> x)
fhn@skeleton@type <- 3L
try(trajectory(fhn))

pompExample(sir)
trajectory(sir,as.data.frame=TRUE) -> x
plot(cases~time,data=x,type='l')

pompExample(gompertz)
gompertz %>% trajectory() -> x
try(gompertz %>%
    pomp(skeleton=map(function(t,x,params,...)c(3,2))) %>%
    trajectory() -> x)
try(gompertz %>% pomp(skeleton=map(function(t,x,params,...)x,delta.t=-1)))

dev.off()
