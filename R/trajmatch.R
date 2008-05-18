traj.match <- function (object, params, est, method = 'Nelder-Mead', gr = NULL, ...) {
  if (!is(object,'pomp'))
    stop("traj.match error: 'object' must be a 'pomp' object")
  if (!is.character(est))
    stop("traj.match error: 'est' must be a character vector containing the names of the parameters to be estimated")
  par.est <- which(names(params)%in%est)
  guess <- params[par.est]
  opt <- optim(
               par=guess,
               fn=function (x) {
                 p <- params
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
  params[par.est] <- opt$par
  opt$params <- params
  opt$par <- NULL
  opt
}
