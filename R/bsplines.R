bspline.basis <- function (x, nbasis, degree = 3, names = NULL) {
  y <- tryCatch(
      .Call(bspline_basis,x,nbasis,degree),
      error = function (e) {
          stop("in ",sQuote("bspline.basis"),": ",conditionMessage(e),call.=FALSE)
      }
  )
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm 
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop("in ",sQuote("bspline.basis"),": length(",sQuote("names"),") must be either 1 or ",nbasis,call.=FALSE)
    }
  }
  y
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1, names = NULL) {
  y <- tryCatch(
      .Call(periodic_bspline_basis,x,nbasis,degree,period),
      error = function (e) {
          stop("in ",sQuote("periodic.bspline.basis"),": ",conditionMessage(e),call.=FALSE)
      }
  )
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm 
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop("in ",sQuote("periodic.bspline.basis"),": length(",
           sQuote("names"),") must be either 1 or ",nbasis,call.=FALSE)
    }
  }
  y
}
