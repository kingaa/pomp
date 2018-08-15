setClassUnion("listies",
  members=c("pompList","abcList","mif2List","pmcmcList","pfilterList"))

setMethod(
  "show",
  signature=signature(object="listies"),
  definition=function (object) {
    y <- as(object,"list")
    names(y) <- names(object)
    show(y)
  }
)

setMethod(
  "print",
  signature=signature(x="listies"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

setMethod(
  "[",
  signature=signature(x="listies"),
  definition=function(x, i, ...) {
    y <- as(x,"list")
    names(y) <- names(x)
    cl <- class(x)
    y <- unlist(y[i])
    if (is.null(y)) {
      list(NULL)
    } else {
      new(cl,y)
    }
  }
)

setMethod(
  "coef",
  signature=signature(object="listies"),
  definition=function(object, ...) {
    do.call(cbind,lapply(object,coef))
  }
)


setMethod(
  "logLik",
  signature=signature(object="listies"),
  definition=function(object, ...) {
    do.call(c,lapply(object,logLik))
  }
)
