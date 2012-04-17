require(pomp.devel)

data(pertussis.sim)

names(pertussis.sim)
x <- lapply(pertussis.sim,as,"data.frame")

print(lapply(x,tail))
