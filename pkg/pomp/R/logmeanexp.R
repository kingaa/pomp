logmeanexp <- function (x) {
  mx <- max(x)
  mx+log(mean(exp(x-mx)))
}
