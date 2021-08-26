##' Ricker model with Poisson observations.
##'
##' \code{ricker} is a \sQuote{pomp} object encoding a stochastic Ricker model
##' with Poisson measurement error.
##'
##' The state process is \eqn{N_{t+1} = r N_{t} \exp(-c N_{t}+e_{t})}{N[t+1] =
##' r N[t] exp(-c N[t]+e[t])}, where the \eqn{e_t}{e[t]} are i.i.d. normal
##' random deviates with zero mean and variance \eqn{\sigma^2}{sigma^2}.  The
##' observed variables \eqn{y_t}{y[t]} are distributed as
##' \eqn{\mathrm{Poisson}(\phi N_t)}{Poisson(phi N[t])}.
##'
##' @return
##' A \sQuote{pomp} object containing the Ricker model and simulated data.
##'
##' @name ricker
##' @docType data
##' @keywords models
##' @family pomp examples
##' @include pomp.R
##' @importFrom utils read.csv2
##'
##' @param r  intrinsic growth rate
##' @param c  density dependence parameter
##' @param sigma environmental process noise s.d.
##' @param phi sampling rate
##' @param N_0 initial condition
##'
##' @example examples/ricker.R
##' @example examples/ricker-bifdiag.R
##'
NULL

##' @rdname ricker
##' @export
ricker <- function (r = exp(3.8), sigma = 0.3, phi = 10, c = 1,
  N_0 = 7)
{

  dat <- '
"time";"y"
0;68
1;2
2;87
3;0
4;12
5;174
6;0
7;0
8;1
9;57
10;11
11;178
12;0
13;1
14;0
15;34
16;72
17;3
18;101
19;0
20;8
21;156
22;0
23;0
24;3
25;93
26;0
27;17
28;121
29;0
30;0
31;19
32;107
33;0
34;4
35;127
36;0
37;1
38;47
39;8
40;117
41;0
42;3
43;82
44;2
45;39
46;70
47;11
48;275
49;0
50;0
'

  pomp(
    data=read.csv2(text=dat),
    times="time", t0=0,
    params=c(r=r,sigma=sigma,phi=phi,c=c,N_0=N_0,e_0=0),
    cfile="ricker_source",
    rprocess=discrete_time(
      step.fun=Csnippet("
      e = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
      N = exp(log(r)+log(N)-c*N+e);"),
      delta.t=1
    ),
    emeasure=Csnippet("
    E_y = phi*N;"),
    rmeasure=Csnippet("
    y = rpois(phi*N);"),
    dmeasure=Csnippet("
    lik = dpois(y,phi*N,give_log);"),
    skeleton=map(Csnippet("
        DN = exp(log(r)+log(N)-c*N);
        De = 0.0;"),
      delta.t=1),
    partrans=parameter_trans(log=c("r","sigma","phi","c","N_0")),
    rinit=Csnippet("
        N = N_0;
        e = 0;"
    ),
    paramnames=c("r","sigma","phi","c","N_0"),
    statenames=c("N","e")
  )

}
