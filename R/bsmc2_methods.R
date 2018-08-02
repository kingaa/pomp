setMethod(
  "eff.sample.size",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

setMethod(
  "logEvidence",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@log.evidence
)

setMethod(
  "cond.logEvidence",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@cond.log.evidence
)

setMethod(
  "prior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@prior[pars,,drop=FALSE]
  })

setMethod(
  "posterior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@post[pars,,drop=FALSE]
  })

setMethod(
  "plot",
  signature(x="bsmcd_pomp"),
  function (x, pars, thin, ...) {
    ep <- paste0("in ",sQuote("plot"),": ")
    if (missing(pars)) pars <- x@est
    pars <- as.character(pars)
    if (length(pars)<1) stop(ep,"no parameters to plot.",call.=FALSE)
    if (missing(thin)) thin <- Inf
    bsmc.plot(
      prior=if (x@transform)
        partrans(x,x@prior,dir="fromEst")
      else
        x@prior,
      post=if (x@transform)
        partrans(x,x@post,dir="fromEst")
      else
        x@post,
      pars=pars,
      thin=thin,
      ...
    )
  }
)

bsmc.plot <- function (prior, post, pars, thin, ...) {
  ep <- paste0("in ",sQuote("plot"),": ")
  p1 <- sample.int(n=ncol(prior),size=min(thin,ncol(prior)))
  p2 <- sample.int(n=ncol(post),size=min(thin,ncol(post)))
  if (!all(pars %in% rownames(prior))) {
    missing <- which(!(pars%in%rownames(prior)))
    stop(ep,"unrecognized parameters: ",
      paste(sQuote(pars[missing]),collapse=","),call.=FALSE)
  }
  prior <- t(prior[pars,,drop=FALSE])
  post <- t(post[pars,,drop=FALSE])
  all <- rbind(prior,post)

  scplot <- function (x, y, ...) { ## prior, posterior pairwise scatterplot
    op <- par(new=TRUE)
    on.exit(par(op))
    i <- which(x[1L]==all[1L,])
    j <- which(y[1L]==all[1L,])
    points(prior[p1,i],prior[p1,j],pch=20,col=rgb(0.85,0.85,0.85,0.1),
      xlim=range(all[,i]),ylim=range(all[,j]))
    points(post[p2,i],post[p2,j],pch=20,col=rgb(0,0,1,0.01))
  }

  dplot <- function (x, ...) { ## marginal posterior histogram
    i <- which(x[1L]==all[1L,])
    d1 <- density(prior[,i])
    d2 <- density(post[,i])
    usr <- par("usr")
    op <- par(usr=c(usr[c(1L,2L)],0,1.5*max(d1$y,d2$y)))
    on.exit(par(op))
    polygon(d1,col=rgb(0.85,0.85,0.85,0.5))
    polygon(d2,col=rgb(0,0,1,0.5))
  }

  if (length(pars) > 1) {
    pairs(all,labels=pars,panel=scplot,diag.panel=dplot)
  }  else {
    d1 <- density(prior[,1])
    d2 <- density(post[,1])
    usr <- par("usr")
    op <- par(usr=c(usr[c(1L,2L)],0,1.5*max(d1$y,d2$y)))
    on.exit(par(op))
    plot(range(all[,1]),range(c(0,d1$y,d2$y)),type='n',
      xlab=pars,ylab="density")
    polygon(d1,col=rgb(0.85,0.85,0.85,0.5))
    polygon(d2,col=rgb(0,0,1,0.5))
  }
}
