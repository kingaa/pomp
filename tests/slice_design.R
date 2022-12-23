png(filename="slice_design-%02d.png",res=100)

library(pomp)
library(dplyr)

## A single 11-point slice through the point c(A=3,B=8,C=0) along the B direction.
x <- slice_design(center=c(A=3,B=8,C=0),B=seq(0,10,by=1))
x |> count(slice) |> as.data.frame()
plot(x)

## Two slices through the same point along the A and C directions.
x <- slice_design(c(A=3,B=8,C=0),A=seq(0,5,by=1),C=seq(0,5,length=11))
x |> count(slice) |> as.data.frame()
plot(x)

try(x <- slice_design(center=c(A=3,C=0),B=seq(0,10,by=1)))
try(x <- slice_design(center=c(A=3),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- slice_design(center=c(3),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- slice_design(center=c(A="3"),B=seq(0,10,by=1),C=c(1,2,3)))
try(x <- slice_design(center=c(A=3,B=2),seq(1,3),A=seq(0,5)))

dev.off()
