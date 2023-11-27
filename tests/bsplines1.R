options(digits=3)
png(filename="bsplines1-%02d.png",res=100)
library(pomp)

x <- seq(0,2,by=0.01)
try(y <- bspline_basis(x,degree=3,nbasis=9,names=c("basis1","basis2")))
y <- bspline_basis(x,degree=3,nbasis=9)
try(y <- bspline_basis(x,degree=3,nbasis=3,names=letters[1:3]))
y <- bspline_basis(x,degree=3,nbasis=12,names=letters[1:12])
y <- bspline_basis(x,degree=3,nbasis=9,names="basis")
y <- bspline_basis(x,degree=3,nbasis=9,names="basis%02d")
matplot(x,y,type='l',ylim=c(0,1.1))
lines(x,apply(y,1,sum),lwd=2)

x <- seq(-1,2,by=0.01)
try(y <- periodic_bspline_basis(x,nbasis=6,names=letters[1:2]))
y <- periodic_bspline_basis(x,nbasis=6,names=tail(letters,6))
y <- periodic_bspline_basis(x,nbasis=5,names="spline")
y <- periodic_bspline_basis(x,nbasis=5,names="spline%d")
matplot(x,y,type='l')

x <- seq(0,1,length=5)
try(bspline_basis(x,degree=-1,nbasis=9))
try(bspline_basis(x,degree=5,nbasis=3))
try(bspline_basis(x,degree=4,nbasis=30,deriv=-5))
try(periodic_bspline_basis(x,degree=-1,nbasis=9))
try(periodic_bspline_basis(x,degree=5,nbasis=4))
try(periodic_bspline_basis(x,degree=4,nbasis=5,deriv=-1))

## now test derivatives
deg <- 5
nb <- 10
dx <- 0.001

x <- seq(0,1,by=dx)
B <- bspline_basis(x,nbasis=nb,degree=deg)
d <- bspline_basis(x,nbasis=nb,degree=deg,deriv=1)
d2 <- bspline_basis(x,nbasis=nb,degree=deg,deriv=2)

B <- apply(B,2,function(x)x-x[1])
dd <- apply(d,2,function(x){y <- diffinv(x); (head(y,-1)+tail(y,-1))/2*dx})
matplot(B,dd,type='l')
abline(a=0,b=1)
stopifnot(all(signif(diag(cor(B,dd)),6)==1))

d <- apply(d,2,function(x) x-x[1])
dd <- apply(d2,2,function(x){y <- diffinv(x); (head(y,-1)+tail(y,-1))/2*dx})
matplot(d,dd,type='l')
abline(a=0,b=1)
stopifnot(all(signif(diag(cor(d,dd)),6)==1))

x <- seq(0,2,by=dx)
B <- periodic_bspline_basis(x,nbasis=nb,degree=deg)
d <- periodic_bspline_basis(x,nbasis=nb,degree=deg,deriv=1)
d2<- periodic_bspline_basis(x,nbasis=nb,degree=deg,deriv=2)

B <- apply(B,2,function(x)x-x[1])
dd <- apply(d,2,function(x){y <- diffinv(x); (head(y,-1)+tail(y,-1))/2*dx})
matplot(B,dd,type='l')
abline(a=0,b=1)
stopifnot(all(signif(diag(cor(B,dd)),6)==1))

d <- apply(d,2,function(x) x-x[1])
dd <- apply(d2,2,function(x){y <- diffinv(x); (head(y,-1)+tail(y,-1))/2*dx})
matplot(d,dd,type='l')
abline(a=0,b=1)
stopifnot(all(signif(diag(cor(d,dd)),6)==1))

B <- bspline_basis(x,degree=0,nbasis=4)
B <- bspline_basis(x,degree=0,nbasis=4,deriv=1)
stopifnot(isTRUE(all(B==0)))
B <- periodic_bspline_basis(x,degree=8,nbasis=30,deriv=11)
stopifnot(isTRUE(all(B==0)))

try(bspline_basis(x,degree=1,nbasis=6,rg=c(4,3)))
try(bspline_basis(x,degree=1,nbasis=6,rg=c(4,4)))
B <- bspline_basis(x,degree=1,nbasis=6,rg=c(-1,3))
matplot(x,B,type='l')

dev.off()
