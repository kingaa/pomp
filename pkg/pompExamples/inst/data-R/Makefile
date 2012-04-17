%.rda: %.R
	Rscript --vanilla $*.R

all: budmoth.sim.rda pertussis.sim.rda
	mv $^ ../../data
