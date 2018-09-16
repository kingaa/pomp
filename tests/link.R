library(pomp2)

cat("double simplefun (double x) { return(x+3); }",file="simplefun.c")
system2(R.home("bin/R"),args=c("CMD","COMPILE","simplefun.c"),
  stdout=NULL,stderr=NULL)

gompertz() -> gompertz

pomp(gompertz,rmeasure=Csnippet("
  double simplefun (double);
  double m = simplefun(X);
  Y = rlnorm(log(m),tau);"),
  statenames="X",paramnames="tau",
  shlib.args="simplefun.o") -> po

x <- simulate(po)

pomp(gompertz,rmeasure=Csnippet("
  double m = simplefun(X);
  Y = rlnorm(log(m),tau);"),
  statenames="X",paramnames="tau",
  shlib.args="simplefun.o",
  globals="double simplefun (double);"
  ) -> po

pomp(gompertz,rmeasure=Csnippet("
  double m = simplefun(X);
  Y = rlnorm(log(m),tau);"),
  statenames="X",paramnames="tau",
  shlib.args="simplefun.o",
  globals=Csnippet("double simplefun (double);")
) -> po

x <- simulate(po)
