.onAttach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  pompExampleDir <- system.file("examples",package="pomp")
  options(pomp.examples=c(exampleDir,pompExampleDir,recursive=TRUE))
}

.onDetach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  pompExampleDir <- system.file("examples",package="pomp")
  exampleDir <- exampleDir[exampleDir!=pompExampleDir]
  options(pomp.examples=exampleDir)
}
