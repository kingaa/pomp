print.pomp <- function (x, ...) {
  print(as(x,'data.frame'))
  invisible(x)
}

setMethod('print','pomp',print.pomp)

show.pomp <- function (object) {
  print.pomp(object)
  invisible(NULL)
}

setMethod('show','pomp',show.pomp)

