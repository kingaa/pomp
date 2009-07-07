traj.match <- function (object, start, est, method = "Nelder-Mead", gr = NULL, ...) {
  if (!is(object,'pomp'))
    stop("traj.match error: ",sQuote("object")," must be a ",sQuote("pomp")," object",call.=FALSE)
  if (is.character(est)) {
    if (!all(est%in%names(start)))
      stop("traj.match error: parameters named in ",sQuote("est")," must exist in ",sQuote("start"),call.=FALSE)
    par.est <- which(names(start)%in%est)
  } else if (is.numeric(est)) {
    est <- as.integer(est)
    if (any((est<1)|(est>length(start))))
      stop("traj.match error: some index in ",sQuote("est")," corresponds to no parameters in ",sQuote("start"),call.=FALSE)
    par.est <- as.integer(est)
  }
  guess <- start[par.est]
  opt <- optim(
               par=guess,
               fn=function (x) {
                 p <- start
                 p[par.est] <- x
                 x <- simulate(object,nsim=1,params=p,states=TRUE)[,,-1,drop=FALSE]
                 d <- dmeasure(
                               object,
                               y=data.array(object),
                               x=x,
                               times=time(object),
                               params=as.matrix(p),
                               log=TRUE
                               )
                 -sum(d)
               },
               gr=gr,
               method=method,
               control=list(...)
               )
  opt$value <- -opt$value
  start[par.est] <- opt$par
  opt$params <- start
  opt$par <- NULL
  opt
}
