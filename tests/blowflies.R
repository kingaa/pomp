library(pomp)

capture.output(
  pompExample(blowflies,envir=NULL) -> flies,
  type="message") -> out
stopifnot(sum(grepl("unrecognized argument",out))==2)

set.seed(599688L)

stopifnot(sum(rinit(flies[[1]]))==11401.5)
x1 <- simulate(flies[[1]])
stopifnot(sum(obs(x1))==499937)
stopifnot(sum(states(x1,"N3"))==496254)
f1 <- pfilter(flies[[1]],Np=1000)
stopifnot(round(logLik(f1),2)==-1466.70)

set.seed(1598688L)
stopifnot(sum(rinit(flies[[2]]))==6037)
x2 <- simulate(flies[[2]])
stopifnot(sum(obs(x2))==561415)
stopifnot(sum(states(x2,"N1"))==560213)
f2 <- pfilter(flies[[2]],Np=1000)
stopifnot(round(logLik(f2),2)==-1477.25)
