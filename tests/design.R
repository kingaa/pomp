library(pomp)
png(filename="design-%02d.png",res=100)

## Sobol' low-discrepancy design
plot(sobolDesign(lower=c(a=0,b=100),upper=c(b=200,a=1),100))

try(sobolDesign(lower=c(a=0,b=100),upper=c(b=200,a=1,q=99),10))
try(sobolDesign(lower=c(0,100),upper=c(b=200,a=1),10))
try(sobolDesign(lower=c(a=0,b=100),upper=c(b=200,c=1),10))

try(sobolDesign(lower=c(a=0,b=100),upper=c(b=200,a=1),2^30+1))
rnames <- sprintf("n%04d",1:5000)
try(sobolDesign(lower=setNames(runif(5000),rnames),
                upper=setNames(runif(5000,min=-1,max=0),rnames),
                100))
x <- sobolDesign(lower=setNames(runif(15),head(rnames,15)),
                 upper=setNames(runif(15,min=1,max=2),head(rnames,15)),
                 100)

## A one-parameter profile design:
x <- profileDesign(p=1:10,lower=c(q=3,a=0,b=0),upper=c(q=5,a=1,b=5),nprof=20)
dim(x)
plot(x)

## A two-parameter profile design:
x <- profileDesign(p=1:10,q=3:5,lower=c(a=0,b=0),upper=c(b=5,a=1),nprof=20)
dim(x)
plot(x)

try(profileDesign(1:10,q=3:5,nprof=10))
try(profileDesign(p=1:10,q=3:5,lower=c(a=0,c=0),upper=c(b=5,a=1),nprof=20))

## A single 11-point slice through the point c(A=3,B=8,C=0) along the B direction.
try(x <- sliceDesign(center=c(A=3,C=0),B=seq(0,10,by=1)))
try(x <- sliceDesign(center=c(A=3),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- sliceDesign(center=c(3),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- sliceDesign(center=c(A="3"),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- sliceDesign(center=c(A=3,B=2),seq(1,3),A=seq(0,5)))
x <- sliceDesign(center=c(A=3,B=8,C=0),B=seq(0,10,by=1))
dim(x)
plot(x)

## Two slices through the same point along the A and C directions.
x <- sliceDesign(c(A=3,B=8,C=0),A=seq(0,5,by=1),C=seq(0,5,length=11))
dim(x)
plot(x)

dev.off()
