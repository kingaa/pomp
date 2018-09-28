verhulst() -> po
plot(po)
plot(simulate(po))
pfilter(po,Np=1000) -> pf
logLik(pf)
spy(po)
