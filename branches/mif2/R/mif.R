## MIF algorithm functions

default.pomp.particles.fun <- function (Np, center, sd, ...) {
	matrix(
			data=rnorm(
					n=Np*length(center),
					mean=center,
					sd=sd
			),
			nrow=length(center),
			ncol=Np,
			dimnames=list(
					names(center),
					NULL
			)
	)
}

mif.cooling <- function (factor, n) {   # default cooling schedule
	alpha <- factor^(n-1)
	list(alpha=alpha,gamma=alpha^2)
}
mif.cooling2 <- function (cooling.fraction, nt, m , n ) {   # cooling schedule for mif2
	#cooling.fraction now is the fraction of cooling on 50 iteration
	cooling.scalar =(50*n*cooling.fraction-1)/(1-cooling.fraction)
	alpha <- (1+cooling.scalar)/(cooling.scalar+nt+n*(m-1))
	list(alpha = alpha, gamma = alpha^2)
}


powerlaw.cooling <- function (init = 1, delta = 0.1, eps = (1-delta)/2, n) {
	m <- init
	if (n <= m) {                         # linear cooling regime
		alpha <- (m-n+1)/m
		gamma <- alpha
	} else {                              # power-law cooling regime
		alpha <- ((n/m)^(delta+eps))/n
		gamma <- ((n/m)^(delta+1))/n/n
	}
	list(alpha=alpha,gamma=gamma)
}

mif.internal <- function (object, Nmif,
		start, pars, ivps,
		particles,
		rw.sd,
		option, cooling.fraction, paramMatrix,
		Np, cooling.factor, var.factor, ic.lag,
		method, tol, max.fail,
		verbose, transform, .ndone) {
	
	transform <- as.logical(transform)
	
	if (length(start)==0)
		stop(
				"mif error: ",sQuote("start")," must be specified if ",
				sQuote("coef(object)")," is NULL",
				call.=FALSE
		)
	
	if (transform)
		start <- partrans(object,start,dir="inverse")
	
	start.names <- names(start)
	if (missing(start.names))
		stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)
	
	if (missing(rw.sd))
		stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
	
	rw.names <- names(rw.sd)
	if (missing(rw.names) || any(rw.sd<0))
		stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
	if (!all(rw.names%in%start.names))
		stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
	rw.names <- names(rw.sd[rw.sd>0])
	if (length(rw.names) == 0)
		stop("mif error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)
	
	if (missing(pars))
		stop("mif error: ",sQuote("pars")," must be specified",call.=FALSE)
	if (missing(ivps))
		stop("mif error: ",sQuote("ivps")," must be specified",call.=FALSE)
	
	if (
			!is.character(pars) ||
			!is.character(ivps) ||
			!all(pars%in%start.names) ||
			!all(ivps%in%start.names) ||
			any(pars%in%ivps) ||
			any(ivps%in%pars) ||
			!all(pars%in%rw.names) ||
			!all(ivps%in%rw.names)
			)
		stop(
				"mif error: ",
				sQuote("pars")," and ",sQuote("ivps"),
				" must be mutually disjoint subsets of ",
				sQuote("names(start)"),
				" and must have a positive random-walk SDs specified in ",
				sQuote("rw.sd"),
				call.=FALSE
		)
	
	if (!all(rw.names%in%c(pars,ivps))) {
		extra.rws <- rw.names[!(rw.names%in%c(pars,ivps))]
		warning(
				"mif warning: the variable(s) ",
				paste(extra.rws,collapse=", "),
				" have positive random-walk SDs specified, but are included in neither ",
				sQuote("pars")," nor ",sQuote("ivps"),
				". These random walk SDs are ignored.",
				call.=FALSE
		)
	}
	rw.sd <- rw.sd[c(pars,ivps)]
	rw.names <- names(rw.sd)
	
	if (missing(particles))
		stop("mif error: ",sQuote("particles")," must be specified",call.=FALSE)
	
	ntimes <- length(time(object))
	if (missing(Np))
		stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
	if (is.function(Np)) {
		Np <- try(
				vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
				silent=FALSE
		)
		if (inherits(Np,"try-error"))
			stop("if ",sQuote("Np")," is a function, it must return a single positive integer")
	}
	if (length(Np)==1)
		Np <- rep(Np,times=ntimes+1)
	else if (length(Np)!=(ntimes+1))
		stop(sQuote("Np")," must have length 1 or length ",ntimes+1)
	if (any(Np<=0))
		stop("number of particles, ",sQuote("Np"),", must always be positive")
	if (!is.numeric(Np))
		stop(sQuote("Np")," must be a number, a vector of numbers, or a function")
	Np <- as.integer(Np)
	
	if (missing(ic.lag))
		stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
	ic.lag <- as.integer(ic.lag)
	if ((length(ic.lag)!=1)||(ic.lag < 1))
		stop("mif error: ",sQuote("ic.lag")," must be a positive integer",call.=FALSE)
	if (ic.lag>ntimes) {
		warning(
				"mif warning: ",sQuote("ic.lag")," = ",ic.lag," > ",ntimes,
				" = length(time(",sQuote("object"),"))",
				" is nonsensical.  Setting ",sQuote("ic.lag")," = ",ntimes,".",
				call.=FALSE
		)
		ic.lag <- length(time(object))
	}
	if ((length(pars)==0)&&(ic.lag<length(time(object)))) {
		warning(
				"mif warning: only IVPs are to be estimated, yet ",sQuote("ic.lag")," = ",ic.lag,
				" < ",ntimes," = length(time(",sQuote("object"),")),",
				" so unnecessary work is to be done.",
				call.=FALSE
		)
	}
	
	if (missing(option) && missing(method))
		stop("mif error: ",sQuote("option")," must be specified",call.=FALSE)
	if (missing(option) && !missing(method) )
	{	option <- method
		warning(sQuote("mif")," warning: ",sQuote("method")," flag is deprecated, use ",sQuote("option"))
	}
	if (missing(cooling.factor)&&(option=="mif2"))	##Default value for the slot cooling.fraction
		cooling.factor=1
	if (missing(cooling.factor)&&(option!="mif2"))
		stop("mif error: ",sQuote("cooling.factor")," must be specified",call.=FALSE)
	if ((option!="mif2")&&((length(cooling.factor)!=1)||(cooling.factor < 0)||(cooling.factor>1)))
		stop("mif error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
	
	if (missing(var.factor))
		stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
	if ((length(var.factor)!=1)||(var.factor < 0))
		stop("mif error: ",sQuote("var.factor")," must be a positive number",call.=FALSE)
	
	if (missing(Nmif))
		stop("mif error: ",sQuote("Nmif")," must be specified",call.=FALSE)
	Nmif <- as.integer(Nmif)
	if (Nmif<0)
		stop("mif error: ",sQuote("Nmif")," must be a positive integer",call.=FALSE)
	if (option=="mif2" && missing(cooling.fraction))
		stop("mif error: ",sQuote("cooling.fraction")," must be specified for method mif2",call.=FALSE)
	if (option=="mif2")
		cooling.fraction <- as.numeric(cooling.fraction)
	if (missing(cooling.fraction)&&(option!="mif2"))	##Default value for the slot cooling.fraction
		cooling.fraction=-1
	
	theta <- start
	
	sigma <- rep(0,length(start))
	names(sigma) <- start.names
	
	rw.sd <- rw.sd[c(pars,ivps)]
	rw.names <- names(rw.sd)
	
	sigma[rw.names] <- rw.sd
	
	conv.rec <- matrix(
			data=NA,
			nrow=Nmif+1,
			ncol=length(theta)+2,
			dimnames=list(
					seq(.ndone,.ndone+Nmif),
					c('loglik','nfail',names(theta))
			)
	)
	conv.rec[1,] <- c(NA,NA,theta)
	
	if (!all(is.finite(theta[c(pars,ivps)]))) {
		stop(
				sQuote("mif"),
				" error: cannot estimate non-finite parameters: ",
				paste(
						c(pars,ivps)[!is.finite(theta[c(pars,ivps)])],
						collapse=","
				),
				call.=FALSE
		)
	}
	
	obj <- as(object,"pomp")
	
	if (Nmif>0)
		tmp.mif <- new("mif",object,particles=particles,Np=Np) # only needed so that we can use the 'particles' method below
	else
		pfp <- obj
	
	for (n in seq_len(Nmif)) {
		##cooling schedule 
		switch(option, mif2={
					cool.sched <- try(mif.cooling2(cooling.fraction, 1 ,  
									.ndone +n,  ntimes), silent = FALSE)
				},
				{ #not mif2	
					cool.sched <- try(mif.cooling(cooling.factor, .ndone + 
											n), silent = FALSE)
				} )
		if (inherits(cool.sched, "try-error")) 
			stop("mif error: cooling schedule error", call. = FALSE)
		sigma.n <- sigma * cool.sched$alpha
		
		P <- try(particles(tmp.mif, Np = Np[1], center = theta, 
						sd = sigma.n * var.factor), silent = FALSE)
		if (inherits(P, "try-error")) 
			stop("mif error: error in ", sQuote("particles"), 
					call. = FALSE)
		## Setting up parameter switch
		switch(option, mif2={
					if(!((n==1)&&(missing(paramMatrix)))) #use paramMatrix if it exists
					{	   P[pars, ] <- paramMatrix[pars,]
					}
					cooling.m=n	
				},
				{
					
				}
			  )
		pfp <- try(pfilter.internal(object = obj, params = P, 
						tol = tol, max.fail = max.fail, pred.mean = (n == 
									Nmif), pred.var = TRUE, filter.mean = TRUE, save.states = FALSE, 
						save.params = FALSE, .rw.sd = sigma.n[pars],cooling.m=cooling.m, cooling.fraction=cooling.fraction, paramMatrix=paramMatrix,option=option,
						verbose = verbose, transform = transform, .ndone=.ndone), silent = FALSE)
		if (inherits(pfp, "try-error")) 
			stop("mif error: error in ", sQuote("pfilter"), 
					call. = FALSE)
		##update switch
		switch(option, mif = {
					v <- pfp$pred.var[pars, , drop = FALSE]
					v1 <- cool.sched$gamma * (1 + var.factor^2) * 
							sigma[pars]^2
					theta.hat <- cbind(theta[pars], pfp$filter.mean[pars, 
									, drop = FALSE])
					theta[pars] <- theta[pars] + colSums(apply(theta.hat, 
									1, diff)/t(v)) * v1
				}, unweighted = {
					theta.hat <- pfp$filter.mean[pars,,drop=FALSE]
					theta[pars] <- rowMeans(theta.hat)
				}, fp = {
					theta.hat <- pfp$filter.mean[pars, ntimes, drop = FALSE]
					theta[pars] <- theta.hat
				}, mif2 ={
					paramMatrix <- pfp$paramMatrix
					theta.hat <- rowMeans(pfp$paramMatrix[pars, , drop = FALSE])
					theta[pars] <- theta.hat
					
				},)
		theta[ivps] <- pfp$filter.mean[ivps, ic.lag]
		conv.rec[n + 1, -c(1, 2)] <- theta
		conv.rec[n, c(1, 2)] <- c(pfp$loglik, pfp$nfail)
		
		if (verbose) 
			cat("MIF iteration ", n, " of ", Nmif, " completed\n")
				
	}
	## back transform the parameter estimate if necessary
	if (transform)
		theta <- partrans(pfp,theta,dir="forward")
	
	new(
			"mif",
			pfp,
			transform=transform,
			params=theta,
			ivps=ivps,
			pars=pars,
			Nmif=Nmif,
			particles=particles,
			var.factor=var.factor,
			ic.lag=ic.lag,
			cooling.factor=cooling.factor,
			random.walk.sd=sigma[rw.names],
			tol=tol,
			conv.rec=conv.rec,
			option=option,
			cooling.fraction = cooling.fraction
	)
}

setGeneric('mif',function(object,...)standardGeneric("mif"))

setMethod(
		"mif",
		signature=signature(object="pomp"),
		function (object, Nmif = 1,
				start,
				pars, ivps = character(0),
				particles, rw.sd,
				Np, ic.lag, var.factor, cooling.factor,
				weighted, option = c("mif","unweighted","fp","mif2"),cooling.fraction,paramMatrix,
				tol = 1e-17, max.fail = 0,
				verbose = getOption("verbose"),
				transform = FALSE, ...) {
			
			transform <- as.logical(transform)
			
			if (missing(start)) start <- coef(object)
			if (missing(rw.sd))
				stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
			if (missing(pars)) {
				rw.names <- names(rw.sd)[rw.sd>0]
				pars <- rw.names[!(rw.names%in%ivps)]
			}
			if (missing(Np))
				stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
			if (missing(ic.lag))
				stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
			if (missing(var.factor))
				stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
			
			if (missing(option)&& missing(method))
				stop("mif error: ",sQuote("option")," must be specified",call.=FALSE)
			if (missing(option) && !missing(method) )
			{	option <- method
				warning(sQuote("mif")," warning: ",sQuote("method")," flag is deprecated, use ",sQuote("option"))
			}
			if (!missing(option))
				option <- match.arg(option)
			if (missing(cooling.factor)&&(option=="mif2"))	##Default value for the slot cooling.fraction
				cooling.factor=1
			if (missing(cooling.factor)&&(option!="mif2"))
				stop("mif error: ",sQuote("cooling.factor")," must be specified",call.=FALSE)
			if ((option!="mif2")&&((length(cooling.factor)!=1)||(cooling.factor < 0)||(cooling.factor>1)))
				stop("mif error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
				
			
			if (option=="mif2" && missing(cooling.fraction))
				stop("mif error: ",sQuote("cooling.fraction")," must be specified for method mif2",call.=FALSE)
			if (!missing(weighted)) {
				warning(sQuote("mif")," warning: ",sQuote("weighted")," flag is deprecated, use ",sQuote("option"))
				if (weighted) {
					if (option!="mif") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "mif"
				} else {
					if (option!="unweighted") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "unweighted"
				}
			}
			
			if (missing(particles)) {         # use default: normal distribution
				particles <- default.pomp.particles.fun
			} else {
				particles <- match.fun(particles)
				if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
					stop(
							"mif error: ",
							sQuote("particles"),
							" must be a function of prototype ",
							sQuote("particles(Np,center,sd,...)"),
							call.=FALSE
					)
			}
			
			mif.internal(
					object=object,
					Nmif=Nmif,
					start=start,
					pars=pars,
					ivps=ivps,
					particles=particles,
					rw.sd=rw.sd,
					Np=Np,
					cooling.factor=cooling.factor,
					var.factor=var.factor,
					ic.lag=ic.lag,
					option=option,
					cooling.fraction = cooling.fraction,
					paramMatrix= paramMatrix,
					tol=tol,
					max.fail=max.fail,
					verbose=verbose,
					transform=transform,
					.ndone=0
			)
			
		}
)


setMethod(
		"mif",
		signature=signature(object="pfilterd.pomp"),
		function (object, Nmif = 1,
				start,
				pars, ivps = character(0),
				particles, rw.sd,
				Np, ic.lag, var.factor, cooling.factor,
				weighted, option = c("mif","unweighted","fp","mif2"),cooling.fraction, paramMatrix,
				tol = 1e-17, max.fail = 0,
				verbose = getOption("verbose"),
				transform = FALSE, ...) {
			
			transform <- as.logical(transform)
			
			if (missing(start)) start <- coef(object)
			if (missing(rw.sd))
				stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
			if (missing(pars)) {
				rw.names <- names(rw.sd)[rw.sd>0]
				pars <- rw.names[!(rw.names%in%ivps)]
			}
			if (missing(Np)) Np <- object@Np
			if (missing(tol)) tol <- object@tol
			if (missing(ic.lag))
				stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
			if (missing(var.factor))
				stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
			
			if (missing(option))
				option <- object@option
			if ((option!="mif2") && missing(cooling.factor))
				cooling.factor <-object@cooling.factor
			
			if (option=="mif2" && missing(cooling.fraction))
				cooling.fraction <- object@cooling.fraction
			if (option=="mif2" && (missing(paramMatrix)))
				paramMatrix <- object@paramMatrix
			
			if (!missing(weighted)) {
				warning(sQuote("mif")," warning: ",sQuote("weighted")," flag is deprecated, use ",sQuote("option"))
				if (weighted) {
					if (option!="mif") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "mif"
				} else {
					if (option!="unweighted") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "unweighted"
				}
			}
			
			if (missing(particles)) {         # use default: normal distribution
				particles <- default.pomp.particles.fun
			} else {
				particles <- match.fun(particles)
				if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
					stop(
							"mif error: ",
							sQuote("particles"),
							" must be a function of prototype ",
							sQuote("particles(Np,center,sd,...)"),
							call.=FALSE
					)
			}
			
			mif.internal(
					object=as(object,"pomp"),
					Nmif=Nmif,
					start=start,
					pars=pars,
					ivps=ivps,
					particles=particles,
					rw.sd=rw.sd,
					Np=Np,
					cooling.factor=cooling.factor,
					var.factor=var.factor,
					ic.lag=ic.lag,
					option=option,
					cooling.fraction=cooling.fraction,
					paramMatrix=paramMatrix,
					tol=tol,
					max.fail=max.fail,
					verbose=verbose,
					transform=transform,
					.ndone=0
			)
		}
)

setMethod(
		"mif",
		signature=signature(object="mif"),
		function (object, Nmif,
				start,
				pars, ivps,
				particles, rw.sd,
				Np, ic.lag, var.factor, cooling.factor,
				weighted, option = c("mif","unweighted","fp","mif2"),cooling.fraction,paramMatrix,
				tol = 1e-17, max.fail = 0,
				verbose = getOption("verbose"),
				transform, ...) {
			
			if (missing(Nmif)) Nmif <- object@Nmif
			if (missing(start)) start <- coef(object)
			if (missing(pars)) pars <- object@pars
			if (missing(ivps)) ivps <- object@ivps
			if (missing(particles)) particles <- object@particles
			if (missing(rw.sd)) rw.sd <- object@random.walk.sd
			if (missing(Np)) Np <- object@Np
			if (missing(ic.lag)) ic.lag <- object@ic.lag
			if (missing(var.factor)) var.factor <- object@var.factor
			if (missing(cooling.factor)) cooling.factor <- object@cooling.factor
			if (missing(tol)) tol <- object@tol
			if (missing(transform)) transform <- object@transform
			transform <- as.logical(transform)
			
			if (missing(option))
				option <- object@option
			if ((option!="mif2") && missing(cooling.factor))
				cooling.factor <-object@cooling.factor
			
			
			if (option=="mif2" && missing(cooling.fraction))
				cooling.fraction <- object@cooling.fraction
			if (option=="mif2" && (missing(paramMatrix)))
				paramMatrix <- object@paramMatrix
			
			if (!missing(weighted)) {
				warning(sQuote("mif")," warning: ",sQuote("weighted")," flag is deprecated, use ",sQuote("option"))
				if (weighted) {
					if (option!="mif") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "mif"
				} else {
					if (option!="unweighted") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "unweighted"
				}
			}
			
			mif.internal(
					object=as(object,"pomp"),
					Nmif=Nmif,
					start=start,
					pars=pars,
					ivps=ivps,
					particles=particles,
					rw.sd=rw.sd,
					Np=Np,
					cooling.factor=cooling.factor,
					var.factor=var.factor,
					ic.lag=ic.lag,
					option=option,
					cooling.fraction=cooling.fraction,
					paramMatrix=paramMatrix,
					tol=tol,
					max.fail=max.fail,
					verbose=verbose,
					transform=transform,
					.ndone=0
			)
		}
)

setMethod(
		'continue',
		signature=signature(object='mif'),
		function (object, Nmif = 1,
				start,
				pars, ivps,
				particles, rw.sd,
				Np, ic.lag, var.factor, cooling.factor,
				weighted, option = c("mif","unweighted","fp","mif2"),cooling.fraction,paramMatrix,
				tol = 1e-17, max.fail = 0,
				verbose = getOption("verbose"),
				transform, ...) {
			
			ndone <- object@Nmif
			if (missing(start)) start <- coef(object)
			if (missing(pars)) pars <- object@pars
			if (missing(ivps)) ivps <- object@ivps
			if (missing(particles)) particles <- object@particles
			if (missing(rw.sd)) rw.sd <- object@random.walk.sd
			if (missing(Np)) Np <- object@Np
			if (missing(ic.lag)) ic.lag <- object@ic.lag
			if (missing(var.factor)) var.factor <- object@var.factor
			if (missing(tol)) tol <- object@tol
			if (missing(transform)) transform <- object@transform
			transform <- as.logical(transform)
			if (missing(option))
				option <- object@option
			if ((option!="mif2") && missing(cooling.factor))
				cooling.factor <-object@cooling.factor
			
			if (option!=object@option && option =="mif2")
				stop("mif error: ",sQuote("option")," for continue should be the same for mif2 option",call.=FALSE)
			if (option=="mif2" && missing(cooling.fraction))
				cooling.fraction <- object@cooling.fraction
			if (option=="mif2" && (missing(paramMatrix)))
				paramMatrix <- object@paramMatrix
			
			if (!missing(weighted)) {
				warning(sQuote("mif")," warning: ",sQuote("weighted")," flag is deprecated, use ",sQuote("option"))
				if (weighted) {
					if (option!="mif") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "mif"
				} else {
					if (option!="unweighted") {
						warning(sQuote("mif")," warning: use of ",sQuote("weighted")," argument overrides choice of ",sQuote("option"))
					}
					option <- "unweighted"
				}
			}
			
			obj <- mif.internal(
					object=as(object,"pomp"),
					Nmif=Nmif,
					start=start,
					pars=pars,
					ivps=ivps,
					particles=particles,
					rw.sd=rw.sd,
					Np=Np,
					cooling.factor=cooling.factor,
					var.factor=var.factor,
					ic.lag=ic.lag,
					option=option,
					cooling.fraction=cooling.fraction,
					paramMatrix=paramMatrix,
					tol=tol,
					max.fail=max.fail,
					verbose=verbose,
					transform=transform,
					.ndone=ndone
			)
			
			object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1,c('loglik','nfail')]
			obj@conv.rec <- rbind(
					object@conv.rec,
					obj@conv.rec[-1,colnames(object@conv.rec)]
			)
			obj@Nmif <- as.integer(ndone+Nmif)
			
			obj
		}
)

mif.profileDesign <- function (object, profile, lower, upper, nprof, ivps, 
		rw.sd, Np, ic.lag, var.factor, cooling.factor,option, cooling.fraction, paramMatrix, ...)
{
	if (missing(profile)) profile <- list()
	if (missing(lower)) lower <- numeric(0)
	if (missing(upper)) upper <- lower
	if (length(lower)!=length(upper))
		stop(sQuote("lower")," and ",sQuote("upper")," must be of the same length")
	pars <- names(lower)
	if (missing(ivps)) ivps <- character(0)
	Np <- as.integer(Np)
	
	pd <- do.call(profileDesign,c(profile,list(lower=lower,upper=upper,nprof=nprof)))
	
	object <- as(object,"pomp")
	
	pp <- coef(object)
	idx <- !(names(pp)%in%names(pd))
	if (any(idx)) pd <- cbind(pd,as.list(pp[idx]))
	
	ans <- vector(mode="list",length=nrow(pd))
	for (k in seq_len(nrow(pd))) {
		ans[[k]] <- list(
				mf=mif(
						object,
						Nmif=0,
						start=unlist(pd[k,]),
						pars=pars,
						ivps=ivps,
						rw.sd=rw.sd,
						Np=Np,
						ic.lag=ic.lag,
						var.factor=var.factor,
						cooling.factor=cooling.factor,
						option=option,
						cooling.fraction=cooling.fraction,
						paramMatrix=paramMatrix,
						...
				)
		)
	}
	
	ans
}
