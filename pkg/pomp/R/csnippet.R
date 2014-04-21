## a class to hold snippets of C code

setClass(
         "Csnippet",
         slots=c(
           text="character"
           ),
         prototype=prototype(
           text=character(0)
           )
         )

Csnippet <- function (text) {
  new(
      "Csnippet",
      text=as.character(text)
      )
}
