library(pomp)

pompExample(gompertz)
set.seed(995868671L)

theta1 <- coef(gompertz)
theta2 <- coef(gompertz,transform=TRUE)
stopifnot(all.equal(partrans(gompertz,theta1,dir="toEstimationScale"),
  theta2))
stopifnot(all.equal(partrans(gompertz,theta2,dir="fromEstimationScale"),
  theta1))
theta3 <- theta2[order(runif(5))]
stopifnot(any(names(theta1)!=names(theta3)))
stopifnot(all.equal(partrans(gompertz,theta3)[names(theta1)],theta1))

theta3 <- theta2[-1]
partrans(gompertz,theta3)
try(partrans(pomp(gompertz,fromEstimationScale=Csnippet("TK = exp(K);"),
  paramnames="K"),theta3))
try(partrans(pomp(gompertz,fromEstimationScale=Csnippet("TK = exp(K);"),
  toEstimationScale=Csnippet("TK = log(K);"),paramnames="K"),theta3))
try(partrans(pomp(gompertz,
  fromEstimationScale=function(params,...)unname(params),
  toEstimationScale=function(params,...)unname(params)),theta3))

pp <- parmat(coef(gompertz),100)
pp["r",] <- runif(100,0,100)
stopifnot(all.equal(partrans(gompertz,partrans(gompertz,pp,dir="to")),pp))
try(partrans(pomp(gompertz,
  fromEstimationScale=function(params,...)unname(params),
  toEstimationScale=function(params,...)unname(params)),pp))
