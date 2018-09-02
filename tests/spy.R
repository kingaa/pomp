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
    paramnames="r",compile=FALSE,cfile="nancy") -> po2
spy(po2)

pompExample(sir)
spy(sir)

pompExample(rw2)
spy(rw2)

pompExample(sir2)
spy(sir2)

try(spy())
try(spy(list()))

