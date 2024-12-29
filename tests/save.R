library(dplyr)
library(tidyr)
library(pomp)
set.seed(1800076828)
ricker() -> po
options(pomp_archive_dir=tempdir())

simulate(po,nsim=20) |>
  coef() |>
  melt() |>
  pivot_wider() |>
  append_data("tmp.csv",overwrite=TRUE)

simulate(po,nsim=20,times=1:3) |>
  as.data.frame() |>
  rename(.id=.L1) |>
  append_data("tmp.csv") -> dat

data.table::fread(file.path(tempdir(),"tmp.csv")) -> dat1

stopifnot(all.equal(dat,dat1))

try(append_data("bob",file="tmp.csv"))
try(append_data("bob",file="tmp.csv",overwrite=TRUE))
