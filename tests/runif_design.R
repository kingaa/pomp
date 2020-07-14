png(filename="runif_design-%02d.png",res=100)
library(pomp)
library(dplyr)
library(magrittr)

set.seed(818859525)

## A 3D random design:
x <- runif_design(lower=c(q=3,a=0,b=0),upper=c(q=5,a=1,b=5),nseq=20)
stopifnot(x %>% count(q) %>% pull(n) %>% unique() %>% equals(1))
plot(x)

## A 1D random design:
x <- runif_design(lower=c(a=1,b=0),upper=c(b=5,a=1),nseq=30)
stopifnot(x %>% count(a) %>% pull(n) %>% unique() %>% equals(30))
plot(x)

try(runif_design(lower=c(),upper=c(),nseq=10))
try(runif_design(lower=c(a=3),upper=c(),nseq=10))
try(runif_design(lower=c(a=3),upper=c(a=1),nseq=10))
try(runif_design(lower=c(a=3),upper=c(b=10),nseq=10))
try(runif_design(lower=c(a=3),upper=c(a=10),nseq=-1))


dev.off()
