require(pomp.devel)

data(budmoth.sim)

names(budmoth.sim)
x <- lapply(budmoth.sim,as,"data.frame")

print(lapply(x,tail))
