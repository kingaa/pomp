pf <- pfilter(gompertz(),Np=1000)	## use 1000 particles

plot(pf)
logLik(pf)
cond.logLik(pf)			## conditional log-likelihoods
eff.sample.size(pf)             ## effective sample size
logLik(pfilter(pf))      	## run it again with 1000 particles

## run it again with 2000 particles
pf <- pfilter(pf,Np=2000,filter.mean=TRUE,filter.traj=TRUE,save.states=TRUE)
fm <- filter.mean(pf)    	## extract the filtering means
ft <- filter.traj(pf)    	## one draw from the smoothing distribution
ss <- saved.states(pf)          ## the latent-state portion of each particle

as(pf,"data.frame") %>% head()
