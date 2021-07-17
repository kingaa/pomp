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
##' an computation is performed and archived in a named file.
##' If the function is called again and the file is present, the computation is not executed.
##' Instead, the results are loaded from the archive.
##' Moreover, via their optional \code{seed} argument, \code{bake} and \code{stew} can control the pseudorandom-number generator (RNG) for greater reproducibility.
##' After the computation is finished, these functions restore the pre-existing RNG state to avoid side effects.
##'
##' The \code{freeze} function doesn't save results, but does set the RNG state to the specified value and restore it after the computation is complete.
##'
##' Both \code{bake} and \code{stew} first test to see whether \code{file} exists.
##' If it does, \code{bake} reads it using \code{\link{readRDS}} and returns the resulting object.
##' By contrast, \code{stew} loads the file using \code{\link{load}} and copies the objects it contains into the user's workspace (or the environment of the call to \code{stew}).
##'
##' If \code{file} does not exist, then both \code{bake} and \code{stew} evaluate the expression \code{expr};
##' they differ in the results that they save.
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
##' @param dependson arbitrary \R object (optional).
##' Variables on which the computation in \code{expr} depends.
##' A hash of these objects will be archived in \code{file}, along with the results of evaluation \code{expr}.
##' When \code{bake} or \code{stew} are called and \code{file} exists, the hash of these objects will be compared against the archived hash;
##' recomputation is forced when these do not match.
##' The dependencies should be specified as unquoted symbols:
##' use a list if there are multiple dependencies.
##' @param info logical.
##' If \code{TRUE}, the \dQuote{ingredients} of the calculation are returned as a list.
##' In the case of \code{bake}, this list is the \dQuote{ingredients} attribute of the returned object.
##' In the case of \code{stew}, this list is a hidden object named \dQuote{.ingredients}, located in the environment within which \code{stew} was called.
##'
##' @inheritParams base::eval
##'
##' @return \code{bake} returns the value of the evaluated expression \code{expr}.
##' Other objects created in the evaluation of \code{expr} are discarded along with the temporary, local environment created for the evaluation.
##'
##' The latter behavior differs from that of \code{stew}, which returns the names of the objects created during the evaluation of \code{expr}.
##' After \code{stew} completes, these objects are copied into the environment in which \code{stew} was called.
##'
##' \code{freeze} returns the value of evaluated expression \code{expr}.
##' However, \code{freeze} evaluates \code{expr} within the parent environment, so other objects created in the evaluation of \code{expr} will therefore exist after \code{freeze} completes.
##'
##' \code{bake} and \code{stew} store information about the time used in evaluating the expression, the state of the random-number generator, and other information, in the archive file.
##' In the case of \code{bake}, this is recorded in the \dQuote{ingredients} attribute (\code{attr(.,"ingredients")});
##' in the \code{stew} case, this is recorded in an object, \dQuote{.ingredients}, in the archive.
##' This information is returned only if \code{info=TRUE}.
##'
##' @author Aaron A. King
##'
##' @example examples/bake.R
##'
NULL

process_dependencies <- function (dependson, envir, ep)
{
  tryCatch(
    digest(eval(dependson,envir=envir)),
    error = function (e) {
      pStop(ep,"cannot compute hash of dependencies: ",
        conditionMessage(e))
    }
  )
}

reload_check <- function (ingredients, code, deps,
  seed, kind, normal.kind, file, ep)
{
  if (is.null(ingredients)) {
    pStop(ep,sQuote(file)," lacks ingredients.")
  }
  identical(code,ingredients$code) &&
    identical(deps,ingredients$dependencies) &&
    identical(seed,ingredients$seed) &&
    identical(kind,ingredients$kind) &&
    identical(normal.kind,ingredients$normal.kind)
}

##' @rdname bake
##' @importFrom digest digest    
##' @export
bake <- function (
  file, expr,
  seed = NULL, kind = NULL, normal.kind = NULL,
  dependson = NULL, info = FALSE
) {
  expr <- substitute(expr)
  code <- digest(deparse(expr))
  deps <- process_dependencies(
    dependson=substitute(dependson),
    envir=parent.frame(),
    ep="bake"
  )
  info <- as.logical(info)
  reload <- file.exists(file)
  if (reload) {
    val <- readRDS(file)
    reload <- reload_check(
      ingredients=attr(val,"ingredients"),
      code=code,deps=deps,
      seed=seed,kind=kind,normal.kind=normal.kind,
      file=file,ep="bake"
    )
    if (!reload) {
      pWarn("bake","recomputing archive ",file,".")
    }
  }
  if (!reload) {
    tmg <- system.time(
      val <- freeze(
        expr,
        seed=seed,
        kind=kind,
        normal.kind=normal.kind,
        envir=parent.frame(1L),
        enclos=parent.frame(2L)
      )
    )
    if (is.null(val)) {
      pWarn("bake","expression evaluates to NULL,",
        " an empty list will be returned.")
      val <- list()
    }
    attr(val,"ingredients") <- list(
      code=code,
      dependencies=deps,
      system.time=tmg,
      seed=seed,
      kind=kind,
      normal.kind=normal.kind
    )
    saveRDS(val,file=file)
  }
  if (!info) {
    attr(val,"ingredients") <- NULL
  }
  val
}

##' @rdname bake
##' @export
stew <- function (
  file, expr,
  seed = NULL, kind = NULL, normal.kind = NULL,
  dependson = NULL, info = FALSE
) {
  expr <- substitute(expr)
  code <- digest(deparse(expr))
  deps <- process_dependencies(
    dependson=substitute(dependson),
    envir=parent.frame(),
    ep="stew"
  )
  info <- as.logical(info)
  reload <- file.exists(file)
  e <- new.env()
  if (reload) {
    objlist <- load(file,envir=e)
    reload <- reload_check(
      ingredients=e$.ingredients,
      code=code,deps=deps,
      seed=seed,kind=kind,normal.kind=normal.kind,
      file=file,ep="stew"
    )
    if (!reload) {
      pWarn("stew","recomputing archive ",file,".")
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
    e$.ingredients <- list(
      code=code,
      dependencies=deps,
      system.time=tmg,
      seed=seed,
      kind=kind,
      normal.kind=normal.kind
    )
    save(list=objects(envir=e,all.names=TRUE),file=file,envir=e)
  }
  objlist <- objects(envir=e,all.names=info)
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
  }
  val
}
