\donttest{ # takes too long for R CMD check
  po <- dacca()
  plot(po)
  ## MLE:
  coef(po)
  plot(simulate(po))
}
