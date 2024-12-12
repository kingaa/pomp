pf <- pfilter(gompertz(),Np=1000)	## use 1000 particles

plot(pf)
logLik(pf)
cond_logLik(pf)			## conditional log-likelihoods
eff_sample_size(pf)             ## effective sample size
logLik(pfilter(pf))      	## run it again with 1000 particles

## run it again with 2000 particles
pf <- pfilter(pf,Np=2000,filter.mean=TRUE,filter.traj=TRUE,save.states="filter")
fm <- filter_mean(pf) ## extract the filtering means
ft <- filter_traj(pf) ## one draw from the smoothing distribution
ss <- saved_states(pf,format="d") ## the latent-state portion of each particle

as(pf,"data.frame") |> head()
