library(pomp)

set.seed(364121008L)
tt <- seq(0,10)
xx <- cbind(x=runif(11),y=runif(11))

.Call("lookup_in_table",ttable=tt,xtable=xx,t=runif(5))
.Call("lookup_in_table",ttable=runif(11),xtable=xx,t=runif(5))
try(.Call("lookup_in_table",ttable=tt[1:3],xtable=xx,t=runif(5)))
