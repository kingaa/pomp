require(pomp)

dat <- read.csv(text='
## Deaths due to plague during an outbreak on the Island of Bombay
## over the period 17 Dec 1905 to 21 July 1906.
##  From Kermack, W. O. & McKendrick, A. G. (1927)
##  A Contribution to the Mathematical Theory of Epidemics
##  Proceedings of the Royal Society of London, Series A, 115: 700--721.
## "date" is date of end of the week (Saturday)
"week","date","deaths"
1,1905-12-23,4
2,1905-12-30,10
3,1906-01-06,15
4,1906-01-13,18
5,1906-01-20,21
6,1906-01-27,31
7,1906-02-03,51
8,1906-02-10,53
9,1906-02-17,97
10,1906-02-24,125
11,1906-03-03,183
12,1906-03-10,292
13,1906-03-17,390
14,1906-03-24,448
15,1906-03-31,641
16,1906-04-07,771
17,1906-04-14,701
18,1906-04-21,696
19,1906-04-28,867
20,1906-05-05,925
21,1906-05-12,801
22,1906-05-19,580
23,1906-05-26,409
24,1906-06-02,351
25,1906-06-09,210
26,1906-06-16,113
27,1906-06-23,65
28,1906-06-30,52
29,1906-07-07,51
30,1906-07-14,39
31,1906-07-21,33
',comment.char="#")

pomp(data=subset(dat,select=c(week,deaths)),
     times="week",
     t0=0,
     params=c(
       beta=2,delta=1.5,y0=0.0004,theta=54,
       sigma=0.02,
       mu=0,gamma=0.2,ratio=10000
     ),
     rprocess=euler.sim(
       step.fun=Csnippet("
                         double X = exp(x);
                         double Y = exp(y);
                         double dx, dy, dn, dW, ito;
                         dx = (mu*(1.0/X-1)+(delta-beta)*Y)*dt;
                         dy = (beta*X+delta*(Y-1)-gamma-mu)*dt;
                         dn = -delta*Y*dt;
                         dW = rnorm(0,sigma*sqrt(dt));
                         ito = 0.5*sigma*sigma*dt;
                         x += dx - beta*Y*(dW-beta*Y*ito);
                         y += dy + beta*X*(dW+beta*X*ito);
                         n += dn;
                         "
       ),
       delta.t=1/24),
     paramnames=c("beta","delta","mu","gamma","sigma","theta","ratio"),
     statenames=c("x","y","n"),
     measurement.model=deaths~nbinom(mu=ratio*exp(y),size=theta),
     logvar=c("beta","delta","ratio","sigma","theta"),
     logitvar=c("y0"),
     parameter.inv.transform=function (params, logvar, logitvar, ...) {
       params[logvar] <- log(params[logvar])
       params[logitvar] <- qlogis(params[logitvar])
       params
     },
     parameter.transform=function (params, logvar, logitvar, ...) {
       params[logvar] <- exp(params[logvar])
       params[logitvar] <- plogis(params[logitvar])
       params
     },
     initializer=function(params, t0, ...) {
       y0 <- unname(params["y0"])
       c(x=log(1-y0),y=log(y0),n=log(1))
     }
     ) -> bbp

c("bbp")
