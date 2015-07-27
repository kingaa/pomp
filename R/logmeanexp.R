logmeanexp <- function (x, se = FALSE) {
  mx <- max(x)
  mean <- mx+log(mean(exp(x-mx)))
  if (se) {
    n <- length(x)
    jk <- vapply(seq_len(n),
                 function(k) logmeanexp(x[-k]),
                 numeric(1))
    c(mean,se=(n-1)*sd(jk)/sqrt(n))
  } else {
    mean
  }
}
