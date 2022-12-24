options(digits=3)

library(pomp)

try(eff_sample_size())
try(eff_sample_size("bob"))

try(filter_mean())
try(filter_mean("bob"))

try(forecast())
try(forecast("bob"))

try(pred_mean())
try(pred_mean("bob"))

try(pred_var())
try(pred_var("bob"))

try(filter_traj())
try(filter_traj("bob"))

try(traces())
try(traces("bob"))

try(continue())
try(continue("bob"))

try(cond_logLik())
try(cond_logLik("bob"))

try(coef())
try(coef("bob"))

try(coef() <- 3)
try(coef("bob") <- 3)

try(logLik())
logLik("bob")
