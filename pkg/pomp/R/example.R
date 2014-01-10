pompExample <- function (example, envir = .GlobalEnv) {
  example <- as.character(substitute(example))
  if (example=="") {
    avlbl <- list.files(
                        path=system.file("examples",package="pomp"),
                        pattern=".+?R$"
                        )
    avlbl <- gsub("\\.R$","",avlbl)
    avlbl
  } else {
    file <- system.file(
                        file.path("examples",paste(example,".R",sep="")),
                        package="pomp"
                        )
    objs <- source(file,local=TRUE)
    for (i in seq_along(objs$value)) {
      assign(objs$value[i],get(objs$value[i]),envir=envir)
    }
    cat("newly created pomp object(s):\n",objs$value,"\n")
    invisible(NULL)
  }
}
