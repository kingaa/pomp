library(pomp)

gompertz() -> gompertz
set.seed(995868671L)

theta1 <- coef(gompertz)
theta2 <- coef(gompertz,transform=TRUE)
stopifnot(all.equal(partrans(gompertz,theta1,dir="toEst"),theta2))
stopifnot(all.equal(partrans(gompertz,theta2,dir="fromEst"),theta1))
theta3 <- theta2[order(runif(5))]
stopifnot(any(names(theta1)!=names(theta3)))
stopifnot(all.equal(partrans(gompertz,theta3)[names(theta1)],theta1))

theta3 <- theta2[-1]
try(partrans(params=theta1))
try(partrans("gompertz",params=theta1))
try(partrans(gompertz,theta3))
try(partrans(pomp(gompertz,partrans=parameter_trans(from=Csnippet("K = exp(T_K);")),
  paramnames="K"),theta3))
try(partrans(pomp(gompertz,partrans=parameter_trans(from=Csnippet("K = exp(T_K);"),
  to=Csnippet("T_K = log(K);")),paramnames="K"),theta3))
try(partrans(pomp(gompertz,partrans=parameter_trans(log="K"),paramnames="K"),
  theta3))
try(partrans(pomp(gompertz,partrans=parameter_trans(from=function(params,...)unname(params),
  to=function(params,...)unname(params))),theta3))

pp <- parmat(coef(gompertz),100)
pp["r",] <- runif(100,0,100)
stopifnot(all.equal(partrans(gompertz,partrans(gompertz,pp,dir="to")),pp))
try(partrans(pomp(gompertz,partrans=parameter_trans(
  from=function(...)unname(list(...)),to=function(...)unname(list(...)))),pp))
try(partrans(pomp(gompertz,partrans=parameter_trans(
  from=function(...)unname(list(...)),to=NULL)),pp))

parameter_trans(toEst=NULL,fromEst=NULL)
try(parameter_trans(toEst=NULL))
try(parameter_trans(fromEst=NULL))

parameter_trans(logit=c("a","b","c"))
pomp:::has_partrans(gompertz)

stopifnot(
  {
    partrans(
      pomp(gompertz,
        partrans=parameter_trans(
          from=function(r,...) c(r=exp(r)),
          to=function(r,...) c(r=log(r))
        )
      ),
      params=parmat(coef(gompertz),2),
      dir="to"
    ) -> p1

    partrans(
      pomp(gompertz,
        partrans=parameter_trans(
          from=function(r,...) c(r=exp(r)),
          to=function(r,...) c(r=log(r))
        )
      ),
      params=p1,
      dir="from"
    ) -> p2

    all.equal(p2,parmat(coef(gompertz),2))

  }
)

ou2() -> ou2
pp <- parmat(coef(ou2),10)
stopifnot(all.equal(partrans(ou2,partrans(ou2,pp,dir="to"),dir="from"),pp))
