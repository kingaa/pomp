library(pomp)

options(verbose=FALSE)

cat("double simplefun (double x) { return(x+3); }",file="simplefun.c")
system2(R.home("bin/R"),args=c("CMD","COMPILE","simplefun.c"))

pompExample(ricker)

pomp(ricker,rmeasure=Csnippet("
  double simplefun (double);
  double m = simplefun(N);
  y = rpois(phi*m);"),
  statenames="N",paramnames="phi",
  shlib.args="simplefun.o") -> po

x <- simulate(po)

pomp(ricker,rmeasure=Csnippet("
  double m = simplefun(N);
  y = rpois(phi*m);"),
  statenames="N",paramnames="phi",
  shlib.args="simplefun.o",
  globals="double simplefun (double);"
  ) -> po

x <- simulate(po)
