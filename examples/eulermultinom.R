## Simulate 5 realizations of Euler-multinomial random variable:

dn <- reulermultinom(5,size=100,rate=c(a=1,b=2,c=3),dt=0.1)
dn

## Compute the probability mass function at each of the 5 realizations:

deulermultinom(x=dn,size=100,rate=c(1,2,3),dt=0.1)

## Compute the expectation of an Euler-multinomial:

eeulermultinom(size=100,rate=c(a=1,b=2,c=3),dt=0.1)

## An Euler-multinomial with overdispersed transitions:

dt <- 0.1
dW <- rgammawn(sigma=0.1,dt=dt)
reulermultinom(5,size=100,rate=c(a=1,b=2,c=3),dt=dW)
