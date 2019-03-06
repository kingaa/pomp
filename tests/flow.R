library(pomp2)

try(flow())
try(flow("bob"))
sir() -> po
try(flow(po))
try(flow(po,t0=0))
try(flow(po,t0=10))
try(flow(po,t0=10,times=1:10))
try(flow(po,t0=10,times=numeric(0)))
try(flow(po,t0=10,times=NULL))
try(flow(po,t0=10,times=30:11))
try(flow(po,t0=10,times=21:30))
try(flow(po,t0=10,times=21:30,params=coef(po)))
flow(po,t0=0,times=1,params=coef(po),x0=rinit(po)) -> x
