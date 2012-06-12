library(pomp)

data(verhulst)

tail(as(verhulst,"data.frame"))
tail(as.data.frame(verhulst))

coef(verhulst,c("n.0","sigma")) <- c(100,0.2)
time(verhulst) <- 1:100

tail(as.data.frame(simulate(verhulst,seed=1066L)))
