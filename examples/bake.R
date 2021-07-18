\dontrun{
  bake(file="example1.rds",{
    x <- runif(1000)
    mean(x)
  })

  bake(file="example1.rds",{
    x <- runif(1000)
    mean(x)
  })

  bake(file="example1.rds",{
    a <- 3
    x <- runif(1000)
    mean(x)
  })

  a <- 5
  b <- 2

  stew(file="example2.rda",
    dependson=list(a,b),{
      x <- runif(10)
      y <- rnorm(n=10,mean=a*x+b,sd=2)
    })

  plot(x,y)

  set.seed(11)
  runif(2)
  freeze(runif(3),seed=5886730)
  runif(2)
  freeze(runif(3),seed=5886730)
  runif(2)

  set.seed(11)
  runif(2)
  runif(2)
  runif(2)

}
