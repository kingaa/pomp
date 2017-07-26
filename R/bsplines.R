bspline.basis <- function (x, nbasis, degree = 3, deriv = 0, names = NULL) {
  ep <- paste0("in ",sQuote("bspline.basis"),": ")
  y <- tryCatch(
      .Call(bspline_basis,x,nbasis,degree,deriv),
      error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
      }
  )
  if (deriv > degree)
    warning(ep,"returning 0 since ",sQuote("deriv")," > ",sQuote("degree"),
      call.=FALSE)
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop(ep,sQuote("names"),"must be of length 1 or ",nbasis,call.=FALSE)
    }
  }
  y
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1,
                                    deriv = 0, names = NULL) {
  ep <- paste0("in ",sQuote("periodic.bspline.basis"),": ")
  y <- tryCatch(
      .Call(periodic_bspline_basis,x,nbasis,degree,period,deriv),
      error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE)
      }
  )
  if (deriv > degree)
    warning(ep,"returning 0 since ",sQuote("deriv")," > ",sQuote("degree"),
      call.=FALSE)
  if (!is.null(names)) {
    if (length(names)==1) {
      nm <- sprintf(names,seq_len(nbasis))
      if (length(unique(nm))!=nbasis)
        nm <- paste(names,seq_len(nbasis),sep=".")
      colnames(y) <- nm
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      stop(ep,sQuote("names"),"must be of length 1 or ",nbasis,call.=FALSE)
    }
  }
  y
}
