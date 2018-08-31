options(digits=3)

library(pomp)

set.seed(770238753)

p <- runif(50)
x <- rcauchy(50,scale=0.1)
stopifnot(
  all.equal(expit(logit(p)),p),
  all.equal(logit(expit(x)),x)
)

Y <- matrix(rcauchy(50,scale=0.1),5,10)
X <- matrix(rexp(50),5,10)
X <- apply(X,2,function(x)x/sum(x))

stopifnot(
  all.equal(apply(apply(X,2,log_barycentric),2,inv_log_barycentric),X),
  all.equal(apply(apply(Y,2,inv_log_barycentric),2,sum),rep(1,ncol(Y))),
  all.equal(apply(apply(apply(apply(Y,2,inv_log_barycentric),2,log_barycentric)-Y,2,range),2,diff),rep(0,ncol(Y)))
)
