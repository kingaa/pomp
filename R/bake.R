##' Bake, stew, and freeze
##'
##' Tools for reproducible computations.
##'
##' @name bake
##' @rdname bake
##' @include package.R
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
##' @param dependson optional.
##' Variables on which the computation in \code{expr} depends.
##' Hashes of these objects will be stored in \code{file}, along with the results of evaluation \code{expr}.
##' When \code{file} exists, hashes of these objects will be compared against the stored values;
##' recomputation is forced when these do not match.
##' These objects should be specified as unquoted symbols, using \code{c} or \code{list} if there are more than one.
##'
##' @inheritParams base::eval
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

##' @rdname bake
##' @importFrom digest digest    
##' @export
bake <- function (
  file, expr,
  seed = NULL, kind = NULL, normal.kind = NULL,
  dependson = NULL
) {
  expr <- substitute(expr)
  reload <- file.exists(file)
  code <- digest(deparse(expr))
  dependson <- as.list(substitute(dependson))
  if (length(dependson)>1) dependson <- dependson[-1L]
  dependson <- lapply(
    dependson,
    all.names,
    functions=TRUE,
    unique=TRUE
  )
  dependson <- unique(c(dependson,recursive=TRUE))
  found <- vapply(
    dependson,
    FUN=exists,
    FUN.VALUE=logical(1),
    envir=parent.frame(),
    inherits=FALSE
  )
  deps <- character(sum(found))
  names(deps) <- sort(dependson[found])
  if (!all(found)) {
    unfound <- dependson[!found]
    pWarn("bake",
      ngettext(length(unfound),"dependency ","dependencies "),
      paste(sQuote(unfound),collapse=","),
      " not found."
    )
  }
  for (n in names(deps)) {
    deps[n] <- eval(
      parse(text=sprintf("digest::digest(%s)",n)),
      envir=parent.frame()
    )
  }
  if (reload) {
    val <- readRDS(file)
    reload <- identical(code,attr(val,"code"))
  }
  if (reload) {
    nn <- names(attr(val,"dependencies"))
    reload <- (length(setdiff(names(deps),nn))==0)
  }
  if (reload) {
    for (n in nn) {
      reload <- reload && identical(deps[n],attr(val,"dependencies")[n])
    }
  }
  if (!reload) {
    tmg <- system.time(
      val <- freeze(
        expr,
        seed=seed,
        kind=kind,
        normal.kind=normal.kind,
        envir=parent.frame(1),
        enclos=parent.frame(2)
      )
    )
    if (is.null(val)) {
      pWarn("bake","expression evaluates to NULL")
      val <- paste0("NULL result returned by ",sQuote("bake"))
    }
    attr(val,"code") <- code
    attr(val,"dependencies") <- deps
    attr(val,"system.time") <- tmg
    saveRDS(val,file=file)
  }
  val
}

##' @rdname bake
##' @export
stew <- function (
  file, expr,
  seed = NULL, kind = NULL, normal.kind = NULL,
  dependson = NULL
) {
  expr <- substitute(expr)
  reload <- file.exists(file)
  code <- digest(deparse(expr))
  dependson <- as.list(substitute(dependson))
  if (length(dependson)>1) dependson <- dependson[-1L]
  dependson <- lapply(
    dependson,
    all.names,
    functions=TRUE,
    unique=TRUE
  )
  dependson <- unique(c(dependson,recursive=TRUE))
  found <- vapply(
    dependson,
    FUN=exists,
    FUN.VALUE=logical(1),
    envir=parent.frame(),
    inherits=FALSE
  )
  deps <- character(sum(found))
  names(deps) <- sort(dependson[found])
  if (!all(found)) {
    unfound <- dependson[!found]
    pWarn("stew",
      ngettext(length(unfound),"dependency ","dependencies "),
      paste(sQuote(unfound),collapse=","),
      " not found."
    )
  }
  for (n in names(deps)) {
    deps[n] <- eval(
      parse(text=sprintf("digest::digest(%s)",n)),
      envir=parent.frame()
    )
  }
  e <- new.env()
  if (reload) {
    objlist <- load(file,envir=e)
    if (is.null(e$.code)) {
      pStop("stew",sQuote(file)," is not a stew.")
    } else {
      reload <- identical(code,e$.code)
    }
  }
  if (reload) {
    if (is.null(e$.dependencies)) {
      pStop("stew",sQuote(file)," is not a stew.")
    } else {
      nn <- names(e$.dependencies)
      reload <- identical(names(deps),nn)
    }
  }
  if (reload) {
    for (n in nn) {
      reload <- reload && identical(deps[n],e$.dependencies[n])
    }
  }
  if (!reload) {
    e <- new.env()
    tmg <- system.time(
      freeze(
        expr,
        envir=e,
        enclos=parent.frame(),
        seed=seed,
        kind=kind,
        normal.kind=normal.kind
      )
    )
    e$.system.time <- tmg
    e$.dependencies <- deps
    e$.code <- code
    ## e$.seed <- attr(val,"seed")
    ## e$.kind <- attr(val,"kind")
    ## e$.normal.kind <- attr(val,"normal.kind")
    save(list=objects(envir=e,all.names=TRUE),file=file,envir=e)
  }
  objlist <- objects(envir=e,all.names=FALSE)
  for (obj in objlist) {
    assign(obj,get(obj,envir=e),envir=parent.frame())
  }
  invisible(objlist)
}

##' @rdname bake
##' @export
freeze <- function (expr, seed = NULL, kind = NULL, normal.kind = NULL,
  envir = parent.frame(),
  enclos =  if(is.list(envir) || is.pairlist(envir))
              parent.frame() else baseenv()
) {
  seed <- as.integer(seed)
  rng.control <- (length(seed) > 0)
  if (rng.control) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed,kind=kind,normal.kind=normal.kind)
  }
  val <- eval(expr,envir=envir,enclos=enclos)
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
