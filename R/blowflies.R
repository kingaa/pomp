##' Nicholson's blowflies.
##'
##' \code{blowflies} is a data frame containing the data from several of Nicholson's classic experiments with the Australian sheep blowfly, \emph{Lucilia cuprina}.
##'
##' \code{blowflies1()} and \code{blowflies2()} construct \sQuote{pomp} objects encoding stochastic delay-difference equation models.
##' The data for these come from "population I", a control culture.
##' The experiment is described on pp. 163--4 of Nicholson (1957).
##' Unlimited quantities of larval food were provided;
##' the adult food supply (ground liver) was constant at 0.4g per day.
##' The data were taken from the table provided by Brillinger et al. (1980).
##'
##' The models are discrete delay equations:
##' \deqn{R(t+1) \sim \mathrm{Poisson}(P N(t-\tau) \exp{(-N(t-\tau)/N_{0})} e(t+1) {\Delta}t)}{R[t+1] ~ Poisson(P N[t-tau] exp(-N[t-tau]/N0) e[t+1] dt)}
##' \deqn{S(t+1) \sim \mathrm{Binomial}(N(t),\exp{(-\delta \epsilon(t+1) {\Delta}t)})}{S[t+1] ~ binomial(N[t],exp(-delta eps[t+1] dt))}
##' \deqn{N(t) = R(t)+S(t)}{N[t]=R[t]+S[t]}
##' where \eqn{e(t)}{e[t]} and \eqn{\epsilon(t)}{eps[t]} are Gamma-distributed i.i.d. random variables
##' with mean 1 and variances \eqn{{\sigma_P^2}/{{\Delta}t}}{sigma.P^2/dt}, \eqn{{\sigma_d^2}/{{\Delta}t}}{sigma.d^2/dt}, respectively.
##' \code{blowflies1} has a timestep (\eqn{{\Delta}t}{dt}) of 1 day; \code{blowflies2} has a timestep of 2 days.
##' The process model in \code{blowflies1} thus corresponds exactly to that studied by Wood (2010).
##' The measurement model in both cases is taken to be
##' \deqn{y(t) \sim  \mathrm{NegBin}(N(t),1/\sigma_y^2)}{y[t] ~ negbin(N[t],1/sigma.y^2),}
##' i.e., the observations are assumed to be negative-binomially distributed with
##' mean \eqn{N(t)}{N[t]} and variance \eqn{N(t)+(\sigma_y N(t))^2}{N[t]+(sigma.y N[t])^2}.
##'
##' Default parameter values are the MLEs as estimated by Ionides (2011).
##'
##' @name blowflies
##' @docType data
##' @family pomp_examples
##' @family pomp_datasets
##' @include pomp.R
##' @importFrom utils read.csv2
##' @importFrom stats approx
##'
##' @return
##' \code{blowflies1} and \code{blowflies2} return \sQuote{pomp} objects containing the actual data and two variants of the model.
##'
##' @references
##'
##' \Nicholson1957
##'
##' \Xia2011
##'
##' \Ionides2011
##'
##' \Wood2010
##'
##' \Gurney1980
##'
##' \Brillinger1980
##'
##' @keywords models pomp_datasets
##' @examples
##'
##' plot(blowflies1())
##' plot(blowflies2())
##'
NULL

## blowfly model, with general dt
## here, set up for dt=1 and dt=2
## dt is hard-coded, and initial values are customized for each dt

## following Xia and Tong, the delay is treated as fixed at 14 days
## Xia and Tong claim to be using tau=8 bidays, but on inspection
## their Euler method is really tau=7 bidays

##' @rdname blowflies
##' @name blowflies1
##'
##' @param P reproduction parameter
##' @param delta death rate
##' @param N0 population scale factor
##' @param sigma.P intensity of \eqn{e} noise
##' @param sigma.d intensity of \eqn{eps} noise
##' @param sigma.y measurement error s.d.
##'
##' @export
blowflies1 <- function (
  P = 3.2838, delta = 0.16073, N0 = 679.94,
  sigma.P = 1.3512, sigma.d = 0.74677, sigma.y = 0.026649
)
{

  dat <- read.csv2(text=blowfly.dat)

  pomp(
    data=dat[dat$day>14,],
    times="day",t0=14,
    params=c(P=P,delta=delta,N0=N0,
      sigma.P=sigma.P,sigma.d=sigma.d,sigma.y=sigma.y),
    cfile="blowfly1_source",
    compile=FALSE,
    rprocess=discrete_time(
      step.fun=Csnippet("
        double *N = &N1;
        int tau = 14, k;
        e = rgammawn(sigma_P,dt)/dt;
        eps = rgammawn(sigma_d,dt)/dt;
        R = rpois(P*N[tau]*exp(-N[tau]/N0)*dt*e);
        S = rbinom(N[0],exp(-delta*dt*eps));
        for (k = tau; k > 0; k--) N[k] = N[k-1];
        N[0] = R+S;"
      ),
      delta.t=1
    ),
    dmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      lik = dnbinom_mu(y,size,N1,give_log);"
    ),
    rmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      y = rnbinom_mu(size,N1);"
    ),
    partrans=parameter_trans(
      log=c("P","delta","N0","sigma.P","sigma.d","sigma.y")
    ),
    paramnames=c("P","N0","delta","sigma.P","sigma.d","sigma.y"),
    statenames=c("N1","R","S","e","eps"),
    y.init=with( ## initial data
      dat,
      approx(x=day,y=y,xout=seq(from=0,to=14,by=1),rule=2)$y
    ),
    ## y.init=c(948, 948, 942, 930, 911, 885, 858, 833.7, 801, 748.3, 676, 589.8, 504, 434.9, 397),
    rinit=function (params, t0, y.init, ...) {
      ntau <- length(y.init)
      n <- y.init[ntau:1]
      names(n) <- paste("N",seq_len(ntau),sep="")
      c(n,R=0,S=0,e=0,eps=0)
    }
  )

}

##' @rdname blowflies
##' @name blowflies2
##' @export
blowflies2 <- function (
  P = 2.7319, delta = 0.17377, N0 = 800.31,
  sigma.P = 1.442, sigma.d = 0.76033, sigma.y = 0.010846)
{

  dat <- read.csv2(text=blowfly.dat)

  pomp(
    data=dat[dat$day>14,],
    times="day",t0=14,
    params=c(P=P,delta=delta,N0=N0,
      sigma.P=sigma.P,sigma.d=sigma.d,sigma.y=sigma.y),
    cfile="blowfly2_source",
    compile=FALSE,
    rprocess=discrete_time(
      step.fun=Csnippet("
        double *N = &N1;
        int tau = 7, k;
        e = rgammawn(sigma_P,dt)/dt;
        eps = rgammawn(sigma_d,dt)/dt;
        R = rpois(P*N[tau]*exp(-N[tau]/N0)*dt*e);
        S = rbinom(N[0],exp(-delta*dt*eps));
        for (k = tau; k > 0; k--) N[k] = N[k-1];
        N[0] = R+S;"
      ),
      delta.t=2
    ),
    dmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      lik = dnbinom_mu(y,size,N1,give_log);"
    ),
    rmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      y = rnbinom_mu(size,N1);"
    ),
    partrans=parameter_trans(
      log=c("P","delta","N0","sigma.P","sigma.d","sigma.y")
    ),
    y.init=with( ## initial data
      dat,
      approx(x=day,y=y,xout=seq(from=0,to=14,by=2),rule=2)$y
    ),
    ## y.init=c(948, 942, 911, 858, 801, 676, 504, 397),
    paramnames=c("P","N0","delta","sigma.P","sigma.d","sigma.y"),
    statenames=c("N1","R","S","e","eps"),
    rinit=function (params, t0, y.init, ...) {
      ntau <- length(y.init)
      n <- y.init[ntau:1]
      names(n) <- paste("N",seq_len(ntau),sep="")
      c(n,R=0,S=0,e=0,eps=0)
    }
  )

}

blowfly.dat <- "
day;y
0;948
2;942
4;911
6;858
8;801
10;676
12;504
14;397
16;248
18;146
20;1801
22;6235
24;5974
26;8921
28;6610
30;5973
32;5673
34;3875
36;2361
38;1352
40;1226
42;912
44;521
46;363
48;229
50;142
52;82
54;542
56;939
58;2431
60;3687
62;4543
64;4535
66;5441
68;4412
70;3022
72;2656
74;1967
76;1295
78;915
80;551
82;313
84;167
86;95
88;93
90;60
92;68
94;5259
96;6673
98;5441
100;3987
102;2952
104;3648
106;4222
108;3889
110;2295
112;1509
114;928
116;739
118;566
120;383
122;274
124;192
126;226
128;519
130;1224
132;2236
134;3818
136;6208
138;5996
140;5789
142;6652
144;7939
146;4868
148;3952
150;2712
152;1734
154;1224
156;703
158;508
160;366
162;279
164;243
166;343
168;761
170;1025
172;1221
174;1600
176;2267
178;3290
180;3471
182;3637
184;3703
186;4876
188;5364
190;4890
192;3029
194;1950
196;1225
198;1076
200;905
202;772
204;628
206;473
208;539
210;825
212;1702
214;2868
216;4473
218;5221
220;6592
222;6400
224;4752
226;3521
228;2719
230;1931
232;1500
234;1082
236;849
238;774
240;864
242;1308
244;1624
246;2224
248;2423
250;2959
252;3547
254;7237
256;5218
258;5311
260;4273
262;3270
264;2281
266;1549
268;1091
270;796
272;610
274;445
276;894
278;1454
280;2262
282;2363
284;3847
286;3876
288;3935
290;3479
292;3415
294;3861
296;3571
298;3113
300;2319
302;1630
304;1297
306;861
308;761
310;659
312;701
314;762
316;1188
318;1778
320;2428
322;3806
324;4519
326;5646
328;4851
330;5374
332;4713
334;7367
336;7236
338;5245
340;3636
342;2417
344;1258
346;766
348;479
350;402
352;248
354;254
356;604
358;1346
360;2342
362;3328
364;3599
366;4081
368;7643
370;7919
372;6098
374;6896
376;5634
378;5134
380;4188
382;3469
384;2442
386;1931
388;1790
390;1722
392;1488
394;1416
396;1369
398;1666
"
