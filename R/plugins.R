setClass(
  "pompPlugin",
  slots=c(
    csnippet='logical',
    slotname='character',
    type='integer',
    step.fn="ANY",
    rate.fn="ANY"
  ),
  prototype=prototype(
    csnippet=FALSE,
    slotname=character(0),
    type=0L,
    step.fn=NULL,
    rate.fn=NULL
  )
)

setClass(
  "onestepRprocessPlugin",
  contains="pompPlugin"
)

setClass(
  "discreteRprocessPlugin",
  contains="pompPlugin",
  slots=c(
    delta.t="numeric"
  )
)

setClass(
  "eulerRprocessPlugin",
  contains="pompPlugin",
  slots=c(
    delta.t="numeric"
  )
)

setClass(
  "gillespieRprocessPlugin",
  contains="pompPlugin",
  slots=c(
    hmax="numeric",
    v="matrix"
  )
)

plugin <- function (object, step.fn, rate.fn) {
  if (missing(object)) {
    new("pompPlugin")
  } else {
    if (!missing(step.fn)) object@step.fn <- step.fn
    if (!missing(rate.fn)) object@rate.fn <- rate.fn
    object
  }
}

onestep.sim <- function (step.fun) {
  new("onestepRprocessPlugin",
      step.fn=step.fun,
      slotname="step.fun",
      csnippet=is(step.fun,"Csnippet"),
      type=1L)
}

discrete.time.sim <- function (step.fun, delta.t = 1) {
  new("discreteRprocessPlugin",
      step.fn=step.fun,
      delta.t=delta.t,
      slotname="step.fun",
      csnippet=is(step.fun,"Csnippet"),
      type=2L)
}

euler.sim <- function (step.fun, delta.t) {
  new("eulerRprocessPlugin",
      step.fn=step.fun,
      delta.t=delta.t,
      slotname="step.fun",
      csnippet=is(step.fun,"Csnippet"),
      type=3L)
}

gillespie.sim <- function (rate.fun, v, d, hmax = Inf) {
  ep <- paste0("in ",sQuote("gillespie.sim")," plugin: ")
  if (!missing(d)) {
    warning("argument ",sQuote("d")," is deprecated; updates to the simulation",
            " algorithm have made it unnecessary", call. = FALSE)
  }
  if (!is.matrix(v)) {
    stop(ep,sQuote("v")," must be a matrix.",
         call.=FALSE)
  }
  if (anyDuplicated(rownames(v))){
    stop(ep,"duplicates in rownames of ",sQuote("v"), call.=FALSE)
  }

  new("gillespieRprocessPlugin",
      rate.fn=rate.fun,
      v=v,
      hmax=hmax,
      slotname="rate.fun",
      csnippet=is(rate.fun,"Csnippet"),
      type=4L)
}

gillespie.hl.sim <- function (..., .pre = "", .post = "", hmax = Inf) {
  ep <- paste0("in ",sQuote("gillespie.hl.sim")," plugin: ")
  args <- list(...)

  for (k in seq_along(args)) {
    if (!is.list(args[[k]]) || length(args[[k]]) != 2) {
      stop(ep,"each of the events should be specified using a length-2 list",
           call.=FALSE)
    }
  }

  codeChunks <- lapply(args, "[[", 1)
  stoich <- lapply(args, "[[", 2)

  checkCode <- function (x) {
    inh <- inherits(x, what = c("Csnippet", "character"), which = TRUE)
    if (!any(inh)) {
      stop(ep,"for each event, the first list-element should be a",
           " C snippet or string.", call.=FALSE)
    }
    if (length(x) != 1){
      stop(ep,"for each event, the length of the first list-element",
           " should be 1.", call.=FALSE)
    }
    as(x,"character")
  }

  codeChunks <- lapply(codeChunks, checkCode)

  tryCatch({
    .pre <- paste(as.character(.pre),collapse="\n")
    .post <- paste(as.character(.post),collapse="\n")
  },
  error = function (e) {
    stop(ep,sQuote(".pre")," and ",sQuote(".post"),
         " must be C snippets or strings.",call.=FALSE)
  })

  for (k in seq_along(stoich)) {
    if (!is.numeric(stoich[[k]]) || is.null(names(stoich[[k]]))) {
      stop(ep,"for each event, the second list-element should be",
           " a named numeric vector", call.=FALSE)
    }
  }

  ## Create C snippet of switch statement
  header <- paste0(.pre, "\nswitch (j) {")
  body <- paste0(
    sprintf("case %d:\n{\n%s\n}\nbreak;\n",seq_along(codeChunks),codeChunks),
    collapse="\n"
  )
  footer <- paste0("default:\nerror(\"unrecognized event %d\",j);\nbreak;\n}\n",.post)
  rate.fn <- Csnippet(paste(header, body, footer, sep="\n"))

  ## Create v matrix
  ## By coercing the vectors to a data frame and then using rbind,
  ## we can ensure that all stoichiometric coefficients for the same
  ## state variables are in the same column even if the vectors in
  ## stoich have differently ordered names. Also, rbind will fail if
  ## the set of variables in each data frame is not the same.
  stoichdf <- lapply(stoich, function (x) data.frame(as.list(x),check.names=FALSE))
  v <- t(data.matrix(do.call(rbind, stoichdf)))
  if (anyDuplicated(rownames(v))){
    stop(ep,"redundant or conflicting stoichiometry.",call.=FALSE)
  }

  new("gillespieRprocessPlugin",
      rate.fn=rate.fn,
      v=v,
      hmax=hmax,
      slotname="rate.fun",
      csnippet=TRUE,
      type=4L)
}

onestep.dens <- function (dens.fun, PACKAGE) {
  warning(sQuote("onestep.dens")," is deprecated and will be removed in a ",
          "forthcoming release.  ","Specify ",sQuote("dprocess")," directly ",
          "using a C snippet or R function.",call.=FALSE)
  if (!missing(PACKAGE))
    warning("in ",sQuote("onestep.dens"),": ",sQuote("PACKAGE")," ignored.",
            call.=FALSE)
  dens.fun
}
