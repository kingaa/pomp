library(pomp)

logmeanexp(rep(10,10),se=TRUE)
logmeanexp(10^seq(8,10),se=TRUE)
logmeanexp(c(1.2e-100,1.4e-100,1.8e-100),se=TRUE)
logmeanexp(c(1e-200,1.5e-200,1.6e-200))
logmeanexp(10^seq(8,10),se=NA)
logmeanexp(10^seq(8,10),se=NULL)
logmeanexp(rep(10,10),ess=TRUE)
logmeanexp(10^seq(8,10),ess=TRUE)
logmeanexp(c(1.2e-100,1.4e-100,1.8e-100),se=TRUE,ess=TRUE)

x1 <- freeze(-rexp(10,rate=0.1),seed=1951379414)
x1 <- x1-max(x1)
stopifnot(
  all.equal(logmeanexp(x1),log(mean(exp(x1))))
)
logmeanexp(x1,se=TRUE)
logmeanexp(x1,ess=TRUE)
logmeanexp(x1,ess=TRUE,se=TRUE)

x2 <- freeze(-rexp(10,rate=10),seed=1951379414)
x2 <- x2-max(x2)
stopifnot(
  all.equal(logmeanexp(x2),log(mean(exp(x2))))
)
logmeanexp(x2,se=TRUE)
logmeanexp(x2,ess=TRUE)
logmeanexp(x2,ess=TRUE,se=TRUE)
