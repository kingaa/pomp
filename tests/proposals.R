options(digits=3)

library(pomp2)
library(magrittr)

pompExample("ou2")
ou2 -> po

set.seed(1178744046L)

try(mvn.diag.rw("bob"))
try(mvn.diag.rw(NULL))
try(mvn.diag.rw())
try(mvn.diag.rw(c(3,2)))
f <- mvn.diag.rw(c(a=3,b=2))
f(c(a=0,b=0))

try(mvn.rw(matrix(0,2,2)))
try(mvn.rw(array(dim=c(2,3),dimnames=list(letters[1:2],LETTERS[1:3]))))
try({m <- diag(3); m[3,3] <- 0; rownames(m) <- colnames(m) <- letters[1:3];
    mvn.rw(m) -> f})
cmat <- matrix(c(1,1,0,1),2,2,dimnames=list(letters[1:2],letters[1:2]))
f <- mvn.rw(cmat)
f(c(a=0,b=0))

try(mvn.rw.adaptive(c(a=1,b=1),scale.start=-2))
try(mvn.rw.adaptive(c(a=1,b=1),scale.cooling=2))
try(mvn.rw.adaptive(c(a=1,b=1),shape.start=-3))
try(mvn.rw.adaptive())
try(mvn.rw.adaptive(rw.sd="bob"))
cmat1 <- matrix(c(1,1,0,1),2,2)
try(mvn.rw.adaptive(rw.var=cmat1))
cmat1 <- matrix(c(1,1,0,1,1,0),2,3,dimnames=list(letters[1:2],letters[1:3]))
try(mvn.rw.adaptive(rw.var=cmat1))
cmat1 <- matrix(c(1,1,0,1),2,2,dimnames=list(letters[1:2],letters[3:4]))
try(mvn.rw.adaptive(rw.var=cmat1))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=-300))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=NA))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=NULL))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=3))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=NA))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=3.2))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=-10))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=3))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=NA))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=Inf))
try(mvn.rw.adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target="bob"))
mvn.rw.adaptive(rw.sd=c(alpha.1=0.1,alpha.3=0.0),
  scale.start=5,scale.cooling=0.1,shape.start=10) -> f1
f1(c(alpha.1=1,alpha.3=1),.n=100,.accepts=1000,verbose=FALSE)

mvn.rw.adaptive(rw.sd=c(alpha.1=0.1,alpha.3=0.1),
  scale.start=5,scale.cooling=0.1,shape.start=10) -> f
options(verbose=TRUE) -> op
capture.output(po %>% pmcmc(Nmcmc=200,Np=100,proposal=f) -> mcmc1) -> out
stopifnot(sum(grepl("proposal covariance matrix:",out))==200)
evalq(covmat.emp,envir=environment(f))
options(op)
