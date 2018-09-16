##' Two-dimensional random-walk process
##'
##' \code{rw2} constructs a \sQuote{pomp} object encoding a 2-D Gaussian random walk.
##'
##' The random-walk process is fully but noisily observed.
##'
##' @name rw2
##' @docType data
##' @keywords datasets models
##' @family pomp examples
##'
##' @return
##' A \sQuote{pomp} object containing simulated data.
##'
##' @examples
##'
##' library(magrittr)
##' library(dplyr)
##' library(tidyr)
##' library(ggplot2)
##'
##' rw2() %>% plot()
##'
##' rw2(s1=1,s2=1,tau=0.1) %>%
##'   simulate(nsim=10,format="d") %>%
##'   ggplot(aes(x=y1,y=y2,group=.id,color=.id))+
##'   geom_path()+
##'   guides(color=FALSE)+
##'   theme_bw()
##'
NULL

##' @rdname rw2
##'
##' @param x1_0,x2_0 initial conditions (i.e., latent state variable values at the zero time \code{t0})
##' @param s1,s2 random walk intensities
##' @param tau observation error s.d.
##' @param t0 zero time
##' @param times observation times
##'
##' @export
rw2 <- function (x1_0 = 0, x2_0 = 0, s1 = 1, s2 = 3, tau = 1,
  times = 1:100, t0 = 0)
{

  simulate(
    times=times,t0=t0,
    params=c(x1_0=x1_0,x2_0=x2_0,s1=s1,s2=s2,tau=tau),
    cfile="rw2_source",
    rprocess = onestep(
      Csnippet("
        x1 = rnorm(x1,s1*sqrt(dt));
        x2 = rnorm(x2,s2*sqrt(dt));"
      )
    ),
    dprocess = Csnippet("
        double sdt = sqrt(t_2 - t_1);
        loglik = dnorm(x1_2,x1_1,s1*sdt,1)+
        dnorm(x2_2,x2_1,s2*sdt,1);"
    ),
    rmeasure=Csnippet("
        y1 = rnorm(x1,tau);
        y2 = rnorm(x2,tau);"
    ),
    dmeasure=Csnippet("
        lik = dnorm(y1,x1,tau,1)+dnorm(y2,x2,tau,1);
        lik = (give_log) ? lik : exp(lik);"
    ),
    statenames=c("x1","x2"),
    obsnames=c("y1","y2"),
    paramnames=c("s1","s2","tau"),
    seed=1376784970L
  )
}
