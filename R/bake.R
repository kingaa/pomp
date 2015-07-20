bake <- function (file, expr) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    val <- eval(expr)
    saveRDS(val,file=file)
    val
  }
}

stew <- function (file, expr) {
  if (file.exists(file)) {
    objlist <- load(file)
    for (obj in objlist)
      assign(obj,get(obj),envir=parent.frame())
  } else {
    expr <- substitute(expr)
    e <- new.env()
    eval(expr,envir=e)
    objlist <- objects(envir=e)
    save(list=objlist,file=file,envir=e)
    for (obj in objlist)
      assign(obj,get(obj,envir=e),envir=parent.frame())
  }
  invisible(objlist)
}
