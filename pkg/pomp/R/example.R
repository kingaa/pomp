pompExample <- function (example) {
  if (missing(example)) {
    avlbl <- list.files(
                        path=system.file("examples",package="pomp"),
                        pattern=".+?R$"
                        )
    avlbl <- gsub("\\.R$","",avlbl)
    cat("available pomp examples:\n",sQuote(avlbl),"\n")
  } else {
    file <- system.file(file.path("examples",paste(example,".R",sep="")),package="pomp")
    objs <- source(file,local=TRUE)
    cat("newly created pomp objects:\n",objs$value,"\n")
  }
  invisible(NULL)
}
