options(digits=3,dplyr.summarise.inform=FALSE)
png(filename="pfilter2-%02d.png",res=100)

library(pomp)
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

set.seed(9994847L)

ou2(tau=5) |>
  window(end=10) |>
  pfilter(
    Np=5000,
    save.states="weighted",
    filter.mean=TRUE,
    pred.mean=TRUE,
    pred.var=TRUE,
    filter.traj=TRUE,
  ) -> pf

pf |> as.data.frame() |> names()
pf |> as.data.frame() |> pivot_longer(-time) |> names()
pf |> forecast() |> melt() |> sapply(class)
pf |> forecast(format="d") |> sapply(class)
pf |> filter_mean() |> melt() |> sapply(class)
pf |> filter_mean(format="d") |> sapply(class)
pf |> pred_mean() |> melt() |> sapply(class)
pf |> pred_mean(format="d") |> sapply(class)
pf |> pred_var() |> melt() |> sapply(class)
pf |> pred_var(format="d") |> sapply(class)
pf |> filter_traj() |> melt() |> sapply(class)
pf |> filter_traj(format="d") |> sapply(class)
pf |> saved_states() |> melt() |> sapply(class)
pf |> saved_states(format="d") |> sapply(class)

c(A=pf,B=pf) -> pfs
pfs |> filter_traj(format="d") |> head()
pfs |> filter_traj(format="d") |> dim()
pfs |> filter_traj(format="a") |> dim()

try(pf |> forecast(format="l"))

pf |> as.data.frame() -> dat0

bind_rows(
  eff.sample.size=data.frame(
    time=seq_along(time(pf)),
    value=eff_sample_size(pf)
  ),
  cond.logLik=data.frame(
    time=seq_along(time(pf)),
    value=cond_logLik(pf)
  ),
  .id="name"
) |>
  pivot_wider() |>
  mutate(
    time=time(pf)[as.integer(time)]
  ) -> dat1
full_join(
  eff_sample_size(pf,format="d"),
  cond_logLik(pf,format="d"),
  by="time"
) -> dat2
stopifnot(
  dat0$eff.sample.size==dat2$eff.sample.size,
  dat0$cond.logLik==dat2$cond.logLik,
  all.equal(dat1,dat2,check.attributes=FALSE)
)

bind_rows(
  forecast=melt(forecast(pf)),
  filter.mean=melt(filter_mean(pf)),
  pred.mean=melt(pred_mean(pf)),
  pred.var=melt(pred_var(pf)),
  .id="type"
) |>
  mutate(
    time=time(pf)[as.integer(time)]
  ) |>
  unite(col=variable,type,variable,sep=".") |>
  pivot_wider(names_from=variable) -> dat1
bind_rows(
  forecast=forecast(pf,format="d"),
  filter.mean=filter_mean(pf,format="d"),
  pred.mean=pred_mean(pf,format="d"),
  pred.var=pred_var(pf,format="d"),
  .id="type"
) |>
  unite(col=variable,type,variable,sep=".") |>
  pivot_wider(names_from=variable) -> dat2
stopifnot(
  all.equal(dat1,dat2,check.attributes=FALSE),
  all.equal(dat0$filter.mean.x1,dat2$filter.mean.x1),
  all.equal(dat0$filter.mean.x2,dat2$filter.mean.x2),
  all.equal(dat0$pred.mean.x1,dat2$pred.mean.x1),
  all.equal(dat0$pred.mean.x2,dat2$pred.mean.x2),
  all.equal(dat0$pred.var.x1,dat2$pred.var.x1),
  all.equal(dat0$pred.var.x2,dat2$pred.var.x2)
)

pf |>
  filter_traj() |>
  melt() |>
  mutate(
    time=time(pf,t0=TRUE)[as.integer(time)]
  ) -> dat1
pf |>
  filter_traj(format="d") -> dat2
stopifnot(
  all.equal(dat1,dat2,check.attributes=FALSE)
)

pf |>
  saved_states(format="l") -> dat1
bind_rows(
  melt(dat1$states),
  melt(dat1$weights)
) |>
  mutate(
    variable=coalesce(variable,".log.weight"),
    time=time(pf)[as.integer(.L1)]
  ) |>
  select(-.L1) |>
  arrange(time,.id) |>
  select(time,.id,variable,value) -> dat1
pf |>
  saved_states(format="d") -> dat2
stopifnot(
  all.equal(dat1,dat2,check.attributes=FALSE)
)

pf |>
  saved_states(format="d") |>
  pivot_wider(names_from="variable") |>
  group_by(time) |>
  summarize(
    p=c(0.05,0.5,0.95),
    x1=wquant(x1,weights=exp(.log.weight),probs=p),
    x2=wquant(x2,weights=exp(.log.weight),probs=p)
  ) |>
  ungroup() |>
  pivot_longer(c(x1,x2)) |>
  pivot_wider(names_from=p) |>
  ggplot(aes(x=time,y=`0.5`,ymin=`0.05`,ymax=`0.95`,group=name))+
  geom_ribbon(alpha=0.5)+
  geom_line(color='red')+
  facet_grid(name~.)+
  labs(y="")

dev.off()
