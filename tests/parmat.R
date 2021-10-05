options(digits=7)
library(pomp)

theta <- c(a=1,b=3,c=4,d=5)
p <- parmat(theta,3)
p
p["b",] <- 1:3
p <- parmat(p,2)
p
theta <- array(
  1:30,dim=c(5,3,2),
  dimnames=list(head(letters,5),head(LETTERS,3),NULL)
)
p <- parmat(theta,2)
p
theta <- array(
  1:30,
  dim=c(5,3,2,1,1,1),
  dimnames=list(head(letters,5),head(LETTERS,3),NULL)
)
q <- parmat(theta,2,names=head(LETTERS,12))
stopifnot(
  all.equal(p,parmat(theta,2)),
  p==q
)
