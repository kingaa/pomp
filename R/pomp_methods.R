## this file contains some basic methods definitions

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
  from="pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(cbind(from@times,t(from@data)))
    names(x) <- c(from@timename,rownames(from@data))
    if (length(from@states)>0) {
      nm <- names(x)
      x <- cbind(x,t(from@states))
      names(x) <- c(nm,rownames(from@states))
    }
    cnames <- get_covariate_names(from@covar)
    if (length(cnames) > 0) {
      nm <- c(names(x),cnames)  # perhaps not strictly necessary (but see issue #56)
      y <- .Call(lookup_in_table,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

as.data.frame.pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

## a simple method to extract the array of states
setMethod(
  "states",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    if (length(object@states)==0) {
      NULL
    } else {
      varnames <- rownames(object@states)
      if (missing(vars)) vars <- varnames
      else if (!all(vars%in%varnames))
        stop("in ",sQuote("states"),": some elements of ",
          sQuote("vars")," correspond to no state variable",call.=FALSE)
      x <- object@states[vars,,drop=FALSE]
      dimnames(x) <- setNames(list(vars,time(object)),
        c("variable",object@timename))
      x
    }
  }
)
