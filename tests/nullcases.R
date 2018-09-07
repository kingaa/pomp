library(pomp2)
library(magrittr)

simulate(times=1:10,t0=0,
         statenames="x",
         obsnames="y",
         params=c(x_0=0),
         rprocess=onestep(function(...)c(x=1)),
         rmeasure=Csnippet("")
) -> s1
states(s1)
obs(s1)

simulate(times=1:10,t0=0,
         statenames="x",
         params=c(x_0=0),
         rprocess=onestep(function(...)c(x=1))
) -> s2
states(s2)
obs(s2)
