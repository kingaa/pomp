options(digits=3)
png(filename="bbs-%02d.png",res=100)

library(pomp)

pompExample(bbs)

plot(bbs)

set.seed(48832734L)

stopifnot(all.equal(coef(bbs),partrans(bbs,coef(bbs,transform=TRUE),dir="from")))
po <- simulate(bbs,seed=48832734L)
plot(po)
pf <- pfilter(bbs,Np=1000)
plot(pf)
tj <- trajectory(bbs)
matplot(t(tj[,1,]),type='l',ylab="")

dev.off()
