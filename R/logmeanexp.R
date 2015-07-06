logmeanexp <- function (x, se = FALSE) {
  mx <- max(x)
  mean <- mx+log(mean(exp(x-mx)))
  if (se) {
    se <- sd(exp(x-mx))/exp(mean-mx)
    c(mean,se=se)
  } else {
    mean
  }
}
