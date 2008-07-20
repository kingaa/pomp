sobol <- function (vars, n) {
  if (!is.list(vars) || is.null(names(vars)))
    stop("sobol error: 'vars' must be a named list")
  if (!all(sapply(vars,function(x)is.numeric(x)&&(length(x)==2))))
    stop("sobol error: each entry in 'vars' must specify a range")
  d <- length(vars)
  x <- .Call("sobol_sequence",as.integer(c(d,n)))
  y <- sapply(
              seq(length=d),
              function (k) {
                vars[[k]][1]+(vars[[k]][2]-vars[[k]][1])*x[k,]
              }
              )
  colnames(y) <- names(vars)
  as.data.frame(y)
}
