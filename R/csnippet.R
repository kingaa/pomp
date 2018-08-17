##' C snippets
##'
##' A class to hold snippets of C code.
##'
##' @name Csnippet
##' @rdname csnippet
##' @aliases as,Csnippet-method Csnippet-class
##' @include pomp-package.R
##' @keywords internal
NULL

setClass(
  "Csnippet",
  slots=c(
    text="character"
  ),
  prototype=prototype(
    text=character(0)
  )
)

##' @name Csnippet
##' @rdname csnippet
Csnippet <- function (text) {
  new("Csnippet",text=as.character(text))
}

##' @name as-csnippet
##' @rdname csnippet
##' @aliases coerce,Csnippet,character-method
setAs(
  from="Csnippet",
  to="character",
  def = function (from) {
    from@text
  }
)
