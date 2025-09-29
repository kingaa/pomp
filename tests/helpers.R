set.seed(901772384)

suppressPackageStartupMessages({
  library(pomp)
  library(tidyr)
})

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

c(
  A=ou2(alpha_1=1),
  B=ou2(tau=2,sigma_1=0)
) |>
  coef(
    pars=c("alpha_1","sigma_3","sigma_1","tau")
  )

c(
  A=ou2(alpha_1=1),
  B=ou2(tau=2,sigma_1=0)
) |>
  coef(
    bob=c("alpha_1","sigma_3","sigma_1","tau")
  )

try(logLik())
logLik("bob")

try(states())
try(states("bob"))
ou2() |>
  states(format="d") |>
  head()
c(A=ou2(),B=gompertz()) |>
  states(format="d") |>
  head()

try(obs())
try(obs("bob"))
ou2() |>
  obs(format="d") |>
  head()
c(A=ou2(),B=gompertz()) |>
  obs(format="d") |>
  head()

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

x <- array(
  1:15,
  dim=c(5,3),
  dimnames=list(X=LETTERS[12:16],letters[1:3])
)
y <- melt(x); y
melt(list(x))
melt(list(list(x)))
melt(list(t(x),x))
names(dimnames(x)) <- c("X","Y")
z1 <- melt(list(a=x,b=x[3:5,]))
z1 <- z1[order(z1$.L1,z1$X,z1$Y),]; z1
z2 <- melt(list(a=x,b=x[3:5,c(2,3,1)]))
z2 <- z2[order(z2$.L1,z2$X,z2$Y),];
z3 <- melt(list(a=t(x),b=x[3:5,c(2,3,1)]))
z3 <- z3[order(z3$.L1,z3$X,z3$Y),names(z2)];
stopifnot(
  all.equal(z1,z2,check.attributes=FALSE),
  all.equal(z1,z3,check.attributes=FALSE)
)

list(
  b=c(a=5,b=2),
  c=array(rnorm(3),dim=3,dimnames=list(name=as.character(1:3))),
  d=array(rnorm(3),dim=3),
  e=array(rnorm(2),dim=2,dimnames=list(name=LETTERS[14:15]))
) |> melt()
