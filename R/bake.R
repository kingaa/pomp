bake <- function (file, expr, seed,
                  kind = NULL, normal.kind = NULL) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    rng.control <- !missing(seed)
    if (rng.control) {
      if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
      save.seed <- get(".Random.seed",envir=.GlobalEnv)
      set.seed(seed,kind=kind,normal.kind=normal.kind)
    }
    val <- eval(expr)
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    saveRDS(val,file=file)
    val
  }
}

stew <- function (file, expr, seed,
                  kind = NULL, normal.kind = NULL) {
  if (file.exists(file)) {
    objlist <- load(file)
    for (obj in objlist)
      assign(obj,get(obj),envir=parent.frame())
  } else {
    rng.control <- !missing(seed)
    expr <- substitute(expr)
    e <- new.env()
    if (rng.control) {
      if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
      save.seed <- get(".Random.seed",envir=.GlobalEnv)
      set.seed(seed,kind=kind,normal.kind=normal.kind)
    }
    eval(expr,envir=e)
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    objlist <- objects(envir=e)
    save(list=objlist,file=file,envir=e)
    for (obj in objlist)
      assign(obj,get(obj,envir=e),envir=parent.frame())
  }
  invisible(objlist)
}

freeze <- function (expr, seed,
                    kind = NULL, normal.kind = NULL) {
  rng.control <- !missing(seed)
  if (rng.control) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed,kind=kind,normal.kind=normal.kind)
  } else warning("seed not set!")
  val <- eval(expr)
  if (rng.control)
    assign(".Random.seed",save.seed,envir=.GlobalEnv)
  val
}
