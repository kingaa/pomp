library(pomp)
set.seed(39596886L)

stopifnot(
  `m1 fail`=rgammawn(n=10000,sigma=1,dt=0.1) %>%
    mean() %>%
    all.equal(0.1,tolerance=0.05),
  `m2 fail`=rgammawn(n=10000,sigma=0.1,dt=0.1) %>%
    mean() %>%
    all.equal(0.1,tolerance=0.05),
  `v1 fail`=rgammawn(n=10000,sigma=1,dt=0.1) %>%
    var() %>%
    all.equal(0.1,tolerance=0.05)
)

