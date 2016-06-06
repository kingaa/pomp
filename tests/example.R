library(pomp)

stopifnot(identical(getOption("pomp.examples"),
                    system.file("examples",package="pomp")))

pompExample()
pompExample(gompertz,show=TRUE)
options(pomp.examples=list(getOption("pomp.examples"),
                           getOption("pomp.examples")))
pompExample()
