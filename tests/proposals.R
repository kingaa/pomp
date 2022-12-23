options(digits=3)

library(pomp)

ou2() -> po

set.seed(1178744046L)

try(mvn_diag_rw("bob"))
try(mvn_diag_rw(NULL))
try(mvn_diag_rw())
try(mvn_diag_rw(c(3,2)))
f <- mvn_diag_rw(c(a=3,b=2))
f(c(a=0,b=0))

try(mvn_rw(matrix(0,2,2)))
try(mvn_rw(array(dim=c(2,3),dimnames=list(letters[1:2],LETTERS[1:3]))))
try({m <- diag(3); m[3,3] <- 0; rownames(m) <- colnames(m) <- letters[1:3];
    mvn_rw(m) -> f})
cmat <- matrix(c(1,1,0,1),2,2,dimnames=list(letters[1:2],letters[1:2]))
f <- mvn_rw(cmat)
f(c(a=0,b=0))

try(mvn_rw_adaptive(c(a=1,b=1),scale.start=-2))
try(mvn_rw_adaptive(c(a=1,b=1),scale.cooling=2))
try(mvn_rw_adaptive(c(a=1,b=1),shape.start=-3))
try(mvn_rw_adaptive())
try(mvn_rw_adaptive(rw.sd="bob"))
cmat1 <- matrix(c(1,1,0,1),2,2)
try(mvn_rw_adaptive(rw.var=cmat1))
cmat1 <- matrix(c(1,1,0,1,1,0),2,3,dimnames=list(letters[1:2],letters[1:3]))
try(mvn_rw_adaptive(rw.var=cmat1))
cmat1 <- matrix(c(1,1,0,1),2,2,dimnames=list(letters[1:2],letters[3:4]))
try(mvn_rw_adaptive(rw.var=cmat1))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=-300))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=NA))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=NULL))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=3))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=NA))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=3.2))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=-10))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=3))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=NA))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target=Inf))
try(mvn_rw_adaptive(rw.var=cmat,scale.start=1000,scale.cooling=0.2,
  shape.start=10,target="bob"))
mvn_rw_adaptive(rw.sd=c(alpha_1=0.1,alpha_3=0.0),
  scale.start=5,scale.cooling=0.1,shape.start=10) -> f1
f1(c(alpha_1=1,alpha_3=1),.n=100,.accepts=1000,verbose=FALSE)

mvn_rw_adaptive(rw.sd=c(alpha_1=0.1,alpha_3=0.1),
  scale.start=5,scale.cooling=0.1,shape.start=10) -> f
options(verbose=TRUE) -> op
capture.output(po |> pmcmc(Nmcmc=200,Np=100,proposal=f) -> mcmc1) -> out
stopifnot(sum(grepl("proposal covariance matrix:",out))==200)
evalq(covmat.emp,envir=environment(f))
options(op)
