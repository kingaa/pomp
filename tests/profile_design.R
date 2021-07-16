png(filename="profile_design-%02d.png",res=100)
library(pomp)
library(dplyr)
library(magrittr)
set.seed(722855899)

## A one-parameter profile design:
x <- profile_design(p=1:10,lower=c(q=3,a=0,b=0),upper=c(q=5,a=1,b=5),nprof=20)
stopifnot(x %>% count(p) %>% pull(n) %>% unique() %>% equals(20))
plot(x)

## A two-parameter profile design:
x <- profile_design(p=1:10,q=3:5,lower=c(a=0,b=0),upper=c(b=5,a=1),nprof=30)
stopifnot(x %>% count(p,q) %>% pull(n) %>% unique() %>% equals(30))
plot(x)

try(profile_design(1:10,q=3:5,nprof=10))
try(profile_design(p=1:10,q=3:5,lower=c(a=0,c=0),upper=c(b=5,a=1),nprof=20))

dev.off()
