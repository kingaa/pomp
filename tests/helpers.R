options(digits=3)

library(pomp)

try(eff.sample.size())
try(eff.sample.size("bob"))

try(filter.mean())
try(filter.mean("bob"))

try(forecast())
try(forecast("bob"))

try(pred.mean())
try(pred.mean("bob"))

try(pred.var())
try(pred.var("bob"))

try(filter.traj())
try(filter.traj("bob"))

try(traces())
try(traces("bob"))

try(continue())
try(continue("bob"))

try(cond.logLik())
try(cond.logLik("bob"))

try(coef())
try(coef("bob"))

try(coef() <- 3)
try(coef("bob") <- 3)

try(logLik())
logLik("bob")
