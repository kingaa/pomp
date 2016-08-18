library(pomp)

png(filename="blowflies-%02d.png",res=100)

set.seed(599688L)

pompExample(blowflies)

init.state(blowflies1)
x1 <- simulate(blowflies1)
plot(x1)
f1 <- pfilter(blowflies1,Np=1000)
logLik(f1)


init.state(blowflies2)
x2 <- simulate(blowflies2)
plot(x2,yax.flip=TRUE)
f2 <- pfilter(blowflies2,Np=1000)
logLik(f2)

dev.off()
