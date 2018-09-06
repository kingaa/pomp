##' @include pomp-package.R

.onAttach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  pompExampleDir <- system.file("examples",package="pomp2")
  options(pomp.examples=c(exampleDir,pompExampleDir,recursive=TRUE))
}

.onDetach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  pompExampleDir <- system.file("examples",package="pomp2")
  exampleDir <- exampleDir[exampleDir!=pompExampleDir]
  options(pomp.examples=exampleDir)
}
