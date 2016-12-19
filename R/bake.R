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
    tmg <- system.time(val <- eval(expr))
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    if (rng.control) {
      attr(val,"seed") <- seed
      attr(val,"kind") <- kind
      attr(val,"normal.kind") <- normal.kind
    }
    attr(val,"system.time") <- tmg
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
    tmg <- system.time(eval(expr,envir=e))
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    objlist <- objects(envir=e)
    save(list=objlist,file=file,envir=e)
    for (obj in objlist)
      assign(obj,get(obj,envir=e),envir=parent.frame())
    if (rng.control) {
      attr(objlist,"seed") <- seed
      attr(objlist,"kind") <- kind
      attr(objlist,"normal.kind") <- normal.kind
    }
    attr(objlist,"system.time") <- tmg
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
  } else warning("in ",sQuote("freeze"),": seed not set!",call.=FALSE)
  tmg <- system.time(val <- eval(expr))
  if (rng.control) {
    assign(".Random.seed",save.seed,envir=.GlobalEnv)
    attr(val,"seed") <- seed
    attr(val,"kind") <- kind
    attr(val,"normal.kind") <- normal.kind
  }
  attr(val,"system.time") <- tmg
  val
}
