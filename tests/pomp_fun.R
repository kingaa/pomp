options(digits=3,keep.source=TRUE)

library(pomp2)

pomp2:::pomp_fun()
pomp2:::pomp_fun(NULL)
try(pomp2:::pomp_fun(36,slotname="silly"))
pomp2:::pomp_fun(function(x,y)x+y,proto=quote(bob(x,y)),slotname="blarp")
try(pomp2:::pomp_fun(function(x)x+y,proto=quote(bob(x,y)),slotname="blarp"))
pomp2:::pomp_fun(function(x)x+y,slotname="blarp") -> f
print(f)
pomp2:::pomp_fun("jehosaphat",PACKAGE="harglebarge",slotname="blurp")
pomp2:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork")
pomp2:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork",Cname="yurp")
pomp2:::pomp_fun(Csnippet("x = x;"),paramnames="x",slotname="blork",Cname="yurp",
  libname="splarm")
pomp2:::pomp_fun(f)

