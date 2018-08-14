options(digits=3,keep.source=TRUE)

library(pomp)
library(magrittr)

pomp:::pomp_fun()
pomp:::pomp_fun(NULL)
try(pomp:::pomp_fun(36,slotname="silly"))
pomp:::pomp_fun(function(x,y)x+y,proto=quote(bob(x,y)),slotname="blarp")
try(pomp:::pomp_fun(function(x)x+y,proto=quote(bob(x,y)),slotname="blarp"))
pomp:::pomp_fun(function(x)x+y,slotname="blarp") -> f
print(f)
pomp:::pomp_fun("jehosaphat",PACKAGE="harglebarge",slotname="blurp")
pomp:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork")
pomp:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork",Cname="yurp")
pomp:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork",Cname="yurp",
  libname="splarm")
pomp:::pomp_fun(f)

