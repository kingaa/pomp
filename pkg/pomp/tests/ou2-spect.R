library(pomp)
pompExample(ou2)

gm1 <- spect.match(ou2,
                  kernel.width=3,
                  detrend="mean",
                  nsim=50,
                  est=c("alpha.1","alpha.4"),
                  method="Nelder-Mead")

gm2 <- spect.match(ou2,
                   kernel.width=3,
                   detrend="mean",
                   nsim=49,
                   est=c("alpha.1","alpha.4"),
                   method="Nelder-Mead")
