\donttest{
  bake(file="example1.rds",{
    x <- runif(1000)
    mean(x)
  })

  stew(file="example2.rda",{
    x <- runif(10)
    y <- rnorm(n=10,mean=3*x+5,sd=2)
  })

  plot(x,y)

  freeze(runif(3),seed=5886730)
  freeze(runif(3),seed=5886730)
}
