options(digits=3)
set.seed(58668844L)

library(pomp)

simulate(times=seq(1,100),t0=0,
  nbasis=9L,
  period=50.0,
  msg="hello!",
  params=setNames(runif(n=9,min=-5,max=5),sprintf("beta%d",1:9)),
  rprocess=euler(
    Csnippet("
      static int first = 1;
      if (first) {
        SEXP Msg = get_userdata(\"msg\");
        char *msg = CHAR(STRING_ELT(Msg,0));
        Rprintf(\"%s\\n\",msg);
        first = 0;
      }
      int nbasis = *(get_userdata_int(\"nbasis\"));
      int degree = 3;
      double period = *(get_userdata_double(\"period\"));
      double dxdt[nbasis];
      periodic_bspline_basis_eval(t,period,degree,nbasis,dxdt);
      x += dt*dot_product(nbasis,dxdt,&beta1);"
    ),delta.t=0.01
  ),
  rmeasure=Csnippet("y = x;"),
  rinit=Csnippet("x = 0;"),
  statenames="x",obsnames="y",paramnames=c("beta1","beta2","beta3")
) -> po

try(po %>%
    simulate(rprocess=onestep(
      Csnippet("
      SEXP Msg = get_userdata(\"bob\");
      char *msg = CHAR(STRING_ELT(Msg,0));
      Rprintf(\"%s\\n\",msg);"))))
try(po %>%
    simulate(rprocess=onestep(
      Csnippet("double nbasis = *(get_userdata_double(\"nbasis\"));"))))
try(po %>%
    simulate(rprocess=onestep(
      Csnippet("double nbasis = *(get_userdata_double(\"bob\"));"))))
try(po %>%
    simulate(rprocess=onestep(
      Csnippet("int nbasis = *(get_userdata_int(\"period\"));"))))
try(po %>%
    simulate(rprocess=onestep(
      Csnippet("int nbasis = *(get_userdata_int(\"bob\"));"))))
try(po %>%
    simulate(rprocess=onestep(
      Csnippet("int nbasis = *(get_userdata_int(\"bob\"));")),
      bob=3))
stopifnot(po %>%
    simulate(rprocess=onestep(
      Csnippet("int nbasis = *(get_userdata_int(\"bob\"));")),
      bob=3L) %>% class() %>% magrittr::equals("pomp"))
