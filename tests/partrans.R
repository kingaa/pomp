library(pomp)

pompExample(gompertz)

pars <- coef(gompertz)

new.trans <- function (params, ...) 
{
  params <- c(params[c("X.0","tau","sigma")], exp(params[c("log.r","log.K")]))
  names(params) <- c("X.0","tau","sigma","r","K")
  params
}

new.inv.trans <- function (params, ...) 
{
  params <- c(params[c("X.0","tau","sigma")], log(params[c("r","K")]))
  names(params) <- c("X.0","tau","sigma","log.r","log.K")
  params
}

po <- pomp(
           gompertz,
           toEstimationScale=new.inv.trans,
           fromEstimationScale=new.trans
           )

coef(po,transform=TRUE) <- new.inv.trans(pars)

stopifnot(identical(new.inv.trans(pars),coef(po,transform=TRUE)))
stopifnot(max(abs(coef(gompertz)-coef(po,names(coef(gompertz)))))<1e-16)

pomp(po,toEstimationScale=function(params,...)runif(length(params)),
     fromEstimationScale=function(params,...)runif(3)) -> po
try(coef(po,transform=TRUE))
try(partrans(po,params=parmat(coef(po),10),dir="to"))
