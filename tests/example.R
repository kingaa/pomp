library(pomp)

stopifnot(identical(getOption("pomp.examples"),
                    system.file("examples",package="pomp")))

pompExample()
pompExample(ricker,show=TRUE)
options(pomp.examples=list(getOption("pomp.examples"),
                           getOption("pomp.examples")))
pompExample()
names(pompExample(ricker,envir=NULL))

try(pompExample(bob))

e <- new.env()
pompExample(ricker,envir=e)
ls(envir=e)

try(pompExample(gompertz,envir="bob"))
