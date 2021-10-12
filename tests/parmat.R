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

try(parmat())
try(parmat("bob"))

theta <- c(a=1,b=3,c=4,d=5)
Theta <- array(theta,dim=length(theta),dimnames=list(names(theta)))
p1 <- parmat(Theta,nrep=2,names=c("A","B"))
Theta <- array(theta,dim=c(length(theta),1),
  dimnames=list(names(theta),NULL))
p2 <- parmat(Theta,nrep=2,names=c("A","B"))
Theta <- array(theta,dim=c(length(theta),1,1,1),
  dimnames=list(names(theta),NULL,NULL,NULL))
p3 <- parmat(Theta,nrep=2,names=c("A","B"))
stopifnot(
  identical(p1,p2),
  identical(p1,p3)
)

try({
  theta <- c(a="tom",b=3,c=4,d=5)
  Theta <- array(theta,dim=length(theta),dimnames=list(names(theta)))
  parmat(Theta,nrep=2,names=c("A","B"))
})

try(
  expand.grid(a=1:3,c="hello",stringsAsFactors=FALSE) %>% parmat()
)
try(
  expand.grid(a=1:3,c="hello",stringsAsFactors=TRUE) %>% parmat()
)
expand.grid(a=1:3,b=1:2) %>% parmat()
expand.grid(a=1:3,b=1:2) %>% parmat(nrep=2)
