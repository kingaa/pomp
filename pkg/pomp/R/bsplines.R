bspline.basis <- function (x, nbasis, degree = 3, names = NULL) {
  y <- .Call(bspline_basis,x,nbasis,degree)
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm 
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop(sQuote("length(names)")," must be either 1 or ",nbasis)
    }
  }
  y
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1, names = NULL) {
  y <- .Call(periodic_bspline_basis,x,nbasis,degree,period)
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm 
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop(sQuote("length(names)")," must be either 1 or ",nbasis)
    }
  }
  y
}
