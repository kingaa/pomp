library(pomp)
set.seed(1178744046L)

pompExample(ou2)

f <- mvn.diag.rw(c(a=10,10))
try(pmcmc(ou2,Nmcmc=2,Np=100,proposal=f))
try(abc(ou2,Nmcmc=2,Np=100,proposal=f,probes=list(probe.mean("y1")),
        scale=1,epsilon=1))
try(mvn.rw(matrix(0,2,2)))
try(mvn.rw(array(dim=c(2,3),dimnames=list(letters[1:2],LETTERS[1:3]))))
try({m <- diag(3); m[3,3] <- 0; rownames(m) <- colnames(m) <- letters[1:3];
    mvn.rw(m)})

try(mvn.rw.adaptive(c(a=1,b=1),scale.start=-2))
try(mvn.rw.adaptive(c(a=1,b=1),scale.cooling=2))
try(mvn.rw.adaptive(c(a=1,b=1),shape.start=-3))
