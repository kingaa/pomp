##' Bake, stew, and freeze
##'
##' Tools for reproducible computations.
##'
##' @name bake
##' @rdname bake
##'
##' @details
##' On cooking shows, recipes requiring lengthy baking or stewing are prepared beforehand.
##' The \code{bake} and \code{stew} functions perform analogously:
##' an computation is performed and stored in a named file.
##' If the function is called again and the file is present, the computation is not executed.
##' Instead, the results are loaded from the file in which they were previously stored.
##' Moreover, via their optional \code{seed} argument, \code{bake} and \code{stew} can control the pseudorandom-number generator (RNG) for greater reproducibility.
##' After the computation is finished, these functions restore the pre-existing RNG state to avoid side effects.
##'
##' The \code{freeze} function doesn't save results, but does set the RNG state to the specified value and restore it after the computation is complete.
##'
##' Both \code{bake} and \code{stew} first test to see whether \code{file} exists.
##' If it does, \code{bake} reads it using \code{\link{readRDS}} and returns the resulting object.
##' By contrast, \code{stew} loads the file using \code{\link{load}} and copies the objects it contains into the user's workspace (or the environment of the call to \code{stew}).
##'
##' If \code{file} does not exist, then both \code{bake} and \code{stew} evaluate the expression \code{expr}; they differ in the results that they save.
##' \code{bake} saves the value of the evaluated expression to \code{file} as a single object.
##' The name of that object is not saved.
##' By contrast, \code{stew} creates a local environment within which \code{expr} is evaluated; all objects in that environment are saved (by name) in \code{file}.
##'
##' @param file Name of the binary data file in which the result will be stored or retrieved, as appropriate.
##' For \code{bake}, this will contain a single object and hence be an RDS file (extension \sQuote{rds});
##' for \code{stew}, this will contain one or more named objects and hence be an RDA file (extension \sQuote{rda}).
##' @param expr Expression to be evaluated.
##' @param seed,kind,normal.kind optional.
##' To set the state and of the RNG.
##' See \code{\link{set.seed}}.
##' The default, \code{seed = NULL}, will not change the RNG state.
##' \code{seed} should be a single integer.
##' See \code{\link{set.seed}}.
##'
##' @return \code{bake} returns the value of the evaluated expression \code{expr}.
##' Other objects created in the evaluation of \code{expr} are discarded along with the temporary, local environment created for the evaluation.
##'
##' The latter behavior differs from that of \code{stew}, which returns the names of the objects created during the evaluation of \code{expr}.
##' After \code{stew} completes, these objects exist in the parent environment (that from which \code{stew} was called).
##'
##' \code{freeze} returns the value of evaluated expression \code{expr}.
##' However, \code{freeze} evaluates \code{expr} within the parent environment, so other objects created in the evaluation of \code{expr} will therefore exist after \code{freeze} completes.
##'
##' \code{bake} and \code{stew} return information about the time used in evaluating the expression.
##' This is recorded in the \code{system.time} attribute of the return value.
##' In addition, if \code{seed} is specified, information about the seed (and the kind of random-number generator used) are stored as attributes of the return value.
##'
##' @author Aaron A. King
##'
##' @example examples/bake.R
##'
NULL

##' @importFrom digest digest    
check_digest <- function (x, dig, code, ...) {
  code2 <- attr(x,"code")
  attr(x,"system.time") <- NULL
  attr(x,"code") <- NULL
  attr(x,"result") <- NULL
  attr(x,"normal.kind") <- NULL
  attr(x,"kind") <- NULL
  attr(x,"seed") <- NULL
  identical(digest(x,...),dig) &&
    identical(code,code2)
}

##' @rdname bake
##' @importFrom digest digest    
##' @export
bake <- function (file, expr, seed = NULL, kind = NULL, normal.kind = NULL) {
  expr <- substitute(expr)
  code <- digest(deparse(expr))
  run <- TRUE
  if (file.exists(file)) {
    val <- readRDS(file)
    run <- !identical(code,attr(val,"code"))
  }
  if (run) {
    seed <- as.integer(seed)
    rng.control <- (length(seed) > 0)
    if (rng.control) {
      if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
      save.seed <- get(".Random.seed",envir=.GlobalEnv)
      set.seed(seed,kind=kind,normal.kind=normal.kind)
    }
    tmg <- system.time(val <- eval(expr))
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    if (is.null(val)) {
      pWarn("bake","expression evaluates to NULL")
      val <- paste0("NULL result returned by ",sQuote("bake"))
    }
    if (rng.control) {
      attr(val,"seed") <- seed
      attr(val,"kind") <- kind
      attr(val,"normal.kind") <- normal.kind
    }
    attr(val,"code") <- code
    attr(val,"result") <- digest(val)
    attr(val,"system.time") <- tmg
    saveRDS(val,file=file)
  }
  val
}

##' @rdname bake
##' @export
stew <- function (file, expr, seed = NULL, kind = NULL, normal.kind = NULL) {
  if (file.exists(file)) {
    objlist <- load(file)
    for (obj in objlist)
      assign(obj,get(obj),envir=parent.frame())
  } else {
    seed <- as.integer(seed)
    rng.control <- (length(seed) > 0)
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

##' @rdname bake
##' @export
freeze <- function (expr, seed = NULL, kind = NULL, normal.kind = NULL) {
  seed <- as.integer(seed)
  rng.control <- (length(seed) > 0)
  if (rng.control) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed,kind=kind,normal.kind=normal.kind)
  }
  val <- eval(expr)
  if (rng.control) {
    assign(".Random.seed",save.seed,envir=.GlobalEnv)
    if (!is.null(val)) {
      attr(val,"seed") <- seed
      attr(val,"kind") <- kind
      attr(val,"normal.kind") <- normal.kind
    }
  }
  val
}
