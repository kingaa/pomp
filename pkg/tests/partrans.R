library(pomp)

data(gompertz)

pars <- coef(gompertz,transform=TRUE)
tpars <- coef(gompertz,transform=FALSE)

new.trans <- function (params, ...) 
{
  params <- c(params[c("X.0","tau","sigma")], log(params[c("r","K")]))
  names(params) <- c("X.0","tau","sigma","log.r","log.K")
  params
}

new.inv.trans <- function (params, ...) 
{
  params <- c(params[c("X.0","tau","sigma")], exp(params[c("log.r","log.K")]))
  names(params) <- c("X.0","tau","sigma","r","K")
  params
}

po <- pomp(gompertz,parameter.transform=new.trans,parameter.inv.transform=new.inv.trans)

stopifnot(identical(pars,coef(po,transform=TRUE,names(pars))))

print(tpars)
print(coef(po))
