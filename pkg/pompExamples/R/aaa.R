.onAttach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  newDir <- system.file("examples",package="pompExamples")
  options(pomp.examples=c(exampleDir,newDir,recursive=TRUE))
}

.onDetach <- function (...) {
  exampleDir <- getOption("pomp.examples")
  newDir <- system.file("examples",package="pompExamples")
  exampleDir <- exampleDir[exampleDir!=newDir]
  options(pomp.examples=exampleDir)
}
