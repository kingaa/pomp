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

parameter_trans <- function (toEst, fromEst) {

  ep <- paste0("in ",sQuote("parameter_trans"),": ")

  ## if one transformation is supplied, then both must be
  c1 <- missing(toEst) || is.null(toEst)
  c2 <- missing(fromEst) || is.null(fromEst)

  if (xor(c1,c2)) {
    stop(ep,"if one of ",sQuote("toEst"),", ",sQuote("fromEst"),
      " is supplied, then so must the other be.",call.=FALSE)
  } else if (c1 || (is(toEst,"pomp_fun") && toEst@mode == -1L)) {
    new("partransPlugin",has=FALSE)
  } else {
    new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
}
