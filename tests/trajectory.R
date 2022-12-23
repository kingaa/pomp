options(digits=3)
png(filename="trajectory-%02d.png",res=100)

library(pomp)
library(ggplot2)

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
y <- trajectory(fhn,params=params,ode_control=list(hmax=0.1),format="a")
invisible(y[,,599:601])
matplot(time(fhn),t(y["V",,]),type='l',lty=1)
plot(y[1,,],y[2,,],type='n')
points(y[1,1,],y[2,1,],pch='.',cex=3,col='black')
points(y[1,2,],y[2,2,],pch='.',cex=3,col='red')

## nonautonomous case
pomp(
  data=fhn,
  times=seq(0,100,by=0.01),t0=0,
  covar=covariate_table(
    i=sin(2*pi*times),
    times=seq(0,101,by=0.1)
  ),
  rinit=function(...) {
    c(V=1,R=0)
  },
  skeleton=vectorfield(
    function(V,R,a,b,c,i,...) {
      c(
        V=c*(V-(V^3)/3-R+i),
        R=(V+a-b*R)/c
      )
    }
  )
) -> fhn1

params <- params[c("a","b","c"),]
invisible(skeleton(fhn1,x,t=c(0,3),params=params))
y <- trajectory(fhn1,params=params,ode_control=list(hmax=0.01),format="a")
y[,,199:201]
matplot(time(fhn1),t(y["V",,]),type='l',lty=1)
plot(y[1,,],y[2,,],type='n')
points(y[1,1,],y[2,1,],pch='.',cex=3,col='black')
points(y[1,2,],y[2,2,],pch='.',cex=3,col='red')

invisible(trajectory(fhn,times=c(1,5),format="a"))
try(trajectory(fhn,times=numeric(0)))
try(trajectory(fhn,times=c(1,1,1)))
try(trajectory(fhn,t0=10))
try(trajectory(fhn,params=c(3,2,1)))
try(trajectory(fhn,params=matrix(c(3,2,1,5),2,2)))
try(trajectory(fhn,params=NULL))
try(trajectory(fhn,params=list(a=3,b=2)))
try(trajectory(fhn,ode_control=list(maxsteps=-1)))
try(trajectory(fhn,ode_control=list(maxsteps=1),verbose=TRUE) -> x)
try(trajectory(pomp(fhn,accumvars="q")) -> x)
fhn@skeleton@type <- 3L
stopifnot(
  {
    trajectory(fhn,format="array") -> x
    sum(is.na(x))==1202
  }
)
try(trajectory("fhn"))
try(trajectory())

sir() -> sir
trajectory(sir,format="data.frame") -> x
plot(cases~time,data=x,type='l')

gompertz() -> gompertz
gompertz |> trajectory(format="a") -> x
gompertz |>
  pomp(
    skeleton=map(function(r,X,Y,K,...){
      c(X=r*X*exp(-X/K),Y=Y+X)
    }),
    accumvars=c("Y"),
    params=c(r=17,X_0=1,Y.0=0,K=100)
  ) -> po3
po3 |>
  trajectory(times=seq(1,1000),format="data.frame") -> dat
plot(X~time,data=dat,subset=(time<100),type='l')
plot(X~Y,data=dat)
gompertz |>
  pomp(accumvars=c("X")) |>
  trajectory(times=seq(1,1000,by=10),format="a") -> x
stopifnot(all(x==0))

try(
  po3 |>
    pomp(skeleton=map(function(...)c(X=1,Y=2,Z=3))) |>
    trajectory(params=c(X_0=1,Y_0=0),format="a")
)

stopifnot(
  po3 |> trajectory(times=seq(0,100,by=5),format="a") |> dim()==c(2,1,21),
  po3 |> pomp(skeleton=NULL) |> trajectory(format="a") |> is.na()
)

trajectory(
  t0=0,times=seq(0,10,by=0.1),
  rinit=function(...)c(x=1,y=0),
  skeleton=vectorfield(function(x,y,t,...)c(x=x,y=200*t)),
  format="d"
) -> dat
plot(x~time,data=dat,col=2,type="l",ylab="")
lines(y~time,data=dat,col=3)

dat1 <- dat[c("time","x","y")]
names(dat1) <- c("time","X","Y")
dat1 |>
  trajectory(
    t0=0,times="time",
    rinit=function(...)c(x=1,y=0),
    skeleton=vectorfield(function(x,y,t,...)c(x=x,y=200*t)),
    format="d"
  ) -> dat2
stopifnot(all.equal(dat,dat2))

ricker() -> po
try(trajectory(po,params=c("A","B")))
p <- parmat(coef(po),3)
colnames(p) <- LETTERS[1:3]
p["r",] <- c(5,10,45)
po |>
  trajectory(params=p,times=1:5,format="array") |>
  dimnames() -> dn
stopifnot(
  names(dn)==c("variable",".id","time"),
  dn$.id==LETTERS[1:3],
  is.null(dn$time)
)
po |>
  trajectory(params=p,times=1:5,format="data.frame") -> x
unique(x$.id)
po |>
  trajectory(format="pomps") |>
  trajectory(params=p,times=1:50,format="pomps") -> pos
states(pos) -> x
obs(pos) -> y
stopifnot(
  names(pos)==LETTERS[1:3],
  length(x)==3,
  length(y)==3
)
pos |>
  as.data.frame() |>
  ggplot(aes(x=time,y=N,group=.id,color=.id))+
  geom_line()+
  facet_wrap(~.id,ncol=1)+
  theme_bw()

dev.off()
