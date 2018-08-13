options(digits=3)

library(pomp)
library(magrittr)

pompExample(gompertz,envir=NULL) %>% extract2(1) -> po

pomp(po,partrans=NULL,bob=3,
  covar=covariate_table(a=0:20,b=0:20,times="a")) -> po1
spy(po1)

po1 %>%
  pomp(partrans=parameter_trans(log="r"),params=NULL,
    rinit=function(params,t0,...)params,
    paramnames="r") -> po2
spy(po2)
