library(pomp)

set.seed(47575684)

examples <- list.files(path=system.file("examples",package="pomp"),pattern=".\\.R$",full.names=TRUE)

for (e in examples) {
  source(e,local=TRUE,echo=TRUE)
}
