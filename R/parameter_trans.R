setClass(
  "partransPlugin",
  slots=c(
    has="logical",
    to="ANY",
    from="ANY"
  ),
  prototype=prototype(
    has=FALSE,
    to=NULL,
    from=NULL
  )
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="NULL",fromEst="NULL"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=FALSE)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="pomp_fun",fromEst="pomp_fun"),
  definition=function(toEst, fromEst, ...) {
    if (toEst@mode == -1L || fromEst@mode == -1L)
      new("partransPlugin",has=FALSE)
    else
      new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="missing",fromEst="missing"),
  definition=function(toEst, fromEst, ..., log, logit, barycentric) {
    if (missing(log) && missing(logit) && missing(barycentric))
      new("partransPlugin",has=FALSE)
    else
      parameter_trans.internal(toEst=NULL,fromEst=NULL,
        log=log,logit=logit,barycentric=barycentric)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="Csnippet",fromEst="Csnippet"),
  definition=function(toEst, fromEst, ..., log, logit, barycentric) {
    if (missing(log) && missing(logit) && missing(barycentric))
      new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
    else
      parameter_trans.internal(toEst=as(toEst,"character"),
        fromEst=as(fromEst,"character"),log=log,logit=logit,
        barycentric=barycentric)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="character",fromEst="character"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="function",fromEst="function"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="ANY",fromEst="missing"),
  definition=function(toEst, fromEst, ...) {
    stop("in ",sQuote("parameter_trans"),": if one of ",sQuote("toEst"),", ",
      sQuote("fromEst")," is supplied, then so must the other be.",call.=FALSE)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="missing",fromEst="ANY"),
  definition=function(toEst, fromEst, ...) {
    stop("in ",sQuote("parameter_trans"),": if one of ",sQuote("toEst"),", ",
      sQuote("fromEst")," is supplied, then so must the other be.",call.=FALSE)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="ANY",fromEst="ANY"),
  definition=function(toEst, fromEst, ...) {
    stop(sQuote("parameter_trans")," not defined for objects of class ",
      sQuote(class(toEst)),", ",sQuote(class(fromEst)),".",call.=FALSE)
  }
)

parameter_trans.internal <- function (toEst, fromEst, ...,
  log, logit, barycentric) {

  #   if (xor(missing(toEst),missing(fromEst)))
  #     stop("in ",sQuote("parameter_trans"),": if one of ",sQuote("toEst"),", ",
  #       sQuote("fromEst")," is supplied, then so must the other be.",call.=FALSE)

  if (missing(toEst)) toEst <- NULL
  if (missing(fromEst)) fromEst <- NULL
  if (missing(log)) log <- NULL
  if (missing(logit)) logit <- NULL
  if (missing(barycentric)) barycentric <- list()

  toEst <- as.character(toEst)
  fromEst <- as.character(fromEst)

  log <- cleanForC(as.character(log))
  logit <- cleanForC(as.character(logit))
  if (is.character(barycentric)) barycentric <- list(barycentric)
  barycentric <- lapply(lapply(barycentric,as.character),cleanForC)

  out1 <- textConnection(object="fromEst",open="a",local=TRUE)
  out2 <- textConnection(object="toEst",open="a",local=TRUE)

  if (length(log) > 0) {
    tpl1 <- "\tT{%v%} = exp({%v%});\n"
    tpl2 <- "\tT{%v%} = log({%v%});\n"
    for (v in log) {
      cat(file=out1,render(tpl1,v=v))
      cat(file=out2,render(tpl2,v=v))
    }
  }

  if (length(logit) > 0) {
    tpl1 <- "\tT{%v%} = expit({%v%});\n"
    tpl2 <- "\tT{%v%} = logit({%v%});\n"
    for (v in logit) {
      cat(file=out1,render(tpl1,v=v))
      cat(file=out2,render(tpl2,v=v))
    }
  }

  if (length(barycentric) > 0) {
    tpl1 <- "\tfrom_log_barycentric(&T{%v%},&{%v%},{%n%});"
    tpl2 <- "\tto_log_barycentric(&T{%v%},&{%v%},{%n%});"
    for (b in barycentric) {
      cat(file=out1,render(tpl1,v=b[1],n=length(b)))
      cat(file=out2,render(tpl2,v=b[1],n=length(b)))
    }
  }

  close(out1)
  close(out2)
  fromEst <- paste(fromEst,collapse="\n")
  toEst <- paste(toEst,collapse="\n")

  new("partransPlugin",has=TRUE,to=Csnippet(toEst),from=Csnippet(fromEst))
}
