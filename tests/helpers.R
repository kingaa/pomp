set.seed(901772384)

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

try(states())
try(states("bob"))

try(obs())
try(obs("bob"))

try(melt())
melt("bob")
x <- data.frame(
  a=letters[1:5],
  b=rnorm(5),
  c=as.integer(1:5),
  d=rexp(5)>0.1
)
try(melt(x))
try(melt(x[1:3]))
try(melt(x[2:4]))
melt(x[2:3])

names(x) <- NULL
melt(x[2:3])

try(
  list(
    a=data.frame(x=1:3,y=rnorm(3)),
    b=1:2,
    c=c("A","B")
  ) |> melt()
)
try(
  list(
    a=data.frame(x=1:3,y=rnorm(3)),
    b=1:2,
    c=c(TRUE,FALSE)
  ) |> melt()
)
