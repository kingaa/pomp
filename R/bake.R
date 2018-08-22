##' Bake, stew, and freeze
##'
##' Tools for reproducible computations.
##'
##' @name bake
##' @rdname bake
##'
##' @details
##' On cooking shows, recipes requiring lengthy baking or stewing are prepared beforehand.  The \code{bake} and \code{stew} functions perform analogously:
##' an computation is performed and stored in a named file.
##' If the function is called again and the file is present, the computation is not executed; rather, the results are loaded from the file in which they were previously stored.
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
##' By contrast, \code{stew} creates a local environment within which \code{expr}is evaluated; all objects in that environment are saved (by name) in \code{file}.
##'
##' @param file Name of the binary data file in which the result will be stored or retrieved, as appropriate.
##' For \code{bake}, this will contain a single object and hence be an RDS file (extension \sQuote{rds});
##' for \code{stew}, this will contain one or more named objects and hence be an RDA file (extension \sQuote{rda}).
##' @param expr Expression to be evaluated.
##' @param seed,kind,normal.kind optional.
##' To set the state and, optionally, kind of RNG used.
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
##' All these functions return information about the time used in evaluating the expression.
##' This is recorded in the \code{system.time} attribute of the return value.
##' In addition, if \code{seed} is specified, information about the seed (and the kind of random-number generator used) are stored as attributes of the return value.
##'
##' @author Aaron A. King
##'
##' @examples
##'
##'   \dontrun{
##'   bake(file="example1.rds",{
##'     x <- runif(1000)
##'     mean(x)
##'   })
##'
##'   stew(file="example2.rda",{
##'     x <- runif(10)
##'     y <- rnorm(n=10,mean=3*x+5,sd=2)
##'   })
##'
##'   plot(x,y)
##'   }
##'
##'   freeze(runif(3),seed=5886730)
##'   freeze(runif(3),seed=5886730)
##'

##' @rdname bake
##' @export
bake <- function (file, expr, seed, kind = NULL, normal.kind = NULL) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    rng.control <- !missing(seed)
    if (missing(seed)) seed <- NULL
    else seed <- as.integer(seed)
    if (length(seed) == 0) seed <- NULL
    if (rng.control) {
      if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
      save.seed <- get(".Random.seed",envir=.GlobalEnv)
      set.seed(seed,kind=kind,normal.kind=normal.kind)
    }
    tmg <- system.time(val <- eval(expr))
    if (rng.control)
      assign(".Random.seed",save.seed,envir=.GlobalEnv)
    if (is.null(val)) {
      warning("in ",sQuote("bake"),": expression evaluates to NULL",call.=FALSE)
      val <- paste0("NULL result returned by ",sQuote("bake"))
    }
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

##' @rdname bake
##' @export
stew <- function (file, expr, seed, kind = NULL, normal.kind = NULL) {
  if (file.exists(file)) {
    objlist <- load(file)
    for (obj in objlist)
      assign(obj,get(obj),envir=parent.frame())
  } else {
    rng.control <- !missing(seed)
    if (missing(seed)) seed <- NULL
    else seed <- as.integer(seed)
    if (length(seed) == 0) seed <- NULL
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
freeze <- function (expr, seed, kind = NULL, normal.kind = NULL) {
  rng.control <- !missing(seed)
  if (missing(seed)) seed <- NULL
  else seed <- as.integer(seed)
  seed <- as.integer(seed)
  if (length(seed) == 0) seed <- NULL
  if (rng.control) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed,kind=kind,normal.kind=normal.kind)
  } else warning("in ",sQuote("freeze"),": seed not set!",call.=FALSE)
  tmg <- system.time(val <- eval(expr))
  if (is.null(val)) {
    warning("in ",sQuote("freeze"),": expression evaluates to NULL",call.=FALSE)
    val <- paste0("NULL result returned by ",sQuote("freeze"))
  }
  if (rng.control) {
    assign(".Random.seed",save.seed,envir=.GlobalEnv)
    attr(val,"seed") <- seed
    attr(val,"kind") <- kind
    attr(val,"normal.kind") <- normal.kind
  }
  attr(val,"system.time") <- tmg
  val
}
