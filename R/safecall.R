## a class to hold unevaluated function calls

setClass(
  "safecall",
  slots=c(
    call="call",
    envir="environment"
  ),
  prototype=prototype(
    call=NULL,
    envir=NULL
  )
)

safecall <- function (...) {
  new("safecall",call=match.call(),envir=parent.frame())
}
