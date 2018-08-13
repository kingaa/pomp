library(pomp)

set.seed(364121008L)
ct <- covariate_table(x=runif(11),y=runif(11),times=seq(0,10))
.Call(pomp:::lookup_in_table,covar=ct,t=runif(5))

try(covariate_table(x=runif(11),y=runif(11),times=seq(10,0)))
try(covariate_table(x=runif(4),y=runif(11),times=seq(0,10)))
try(covariate_table(x=runif(11),y=runif(11),times=seq(0,3)))
covariate_table()
try(covariate_table(times=1:10))
try(covariate_table(a=1:10,times=1))
try(covariate_table(a=1:10,times="a"))
