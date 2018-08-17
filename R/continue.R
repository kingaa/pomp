##' Continue an iterative calculation
##'
##' Continue an iterative computation where it left off.
##'
##' @name continue
##' @include mif2.R abc.R pmcmc.R
##' @aliases continue continue,missing-method continue,ANY-method
##' @rdname continue
##'
##' @param object the result of an iterative \pkg{pomp} computation
##' @param \dots additional arguments will be passed to the underlying method.
##' This allows one to modify parameters used in the original computations.
NULL

##' @name continue
##' @rdname continue
setGeneric(
  "continue",
  function (object, ...)
    standardGeneric("continue")
)

setMethod(
  "continue",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("continue"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "continue",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("continue")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name continue-abcd_pomp
##' @aliases continue,abcd_pomp-method
##' @rdname continue
##'
##' @param Nabc positive integer; number of additional iterations to perform
##'
setMethod(
  "continue",
  signature=signature(object="abcd_pomp"),
  definition=function (object, Nabc = 1, ...) {

    ndone <- object@Nabc
    accepts <- object@accepts

    obj <- abc(
      object=object,
      Nabc=Nabc,
      .ndone=ndone,
      .accepts=accepts,
      ...
    )

    obj@traces <- rbind(
      object@traces[,colnames(obj@traces)],
      obj@traces[-1,]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    obj@Nabc <- as.integer(ndone+Nabc)
    obj@accepts <- as.integer(accepts+obj@accepts)

    obj
  }
)

##' @name continue-mif2d_pomp
##' @aliases continue,mif2d_pomp-method
##' @rdname continue
##'
##' @param Nmif positive integer; number of additional iterations to perform
##'
setMethod(
  "continue",
  signature=signature(object="mif2d_pomp"),
  definition = function (object, Nmif = 1, ...) {

    ndone <- object@Nmif

    obj <- mif2(object,Nmif=Nmif,...,
      .ndone=ndone,.paramMatrix=object@paramMatrix)

    object@traces[ndone+1,c('loglik','nfail')] <- obj@traces[1L,c('loglik','nfail')]
    obj@traces <- rbind(
      object@traces,
      obj@traces[-1L,colnames(object@traces)]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    obj@Nmif <- as.integer(ndone+Nmif)

    obj
  }
)

##' @name continue-pmcmcd_pomp
##' @aliases continue,pmcmcd_pomp-method
##' @rdname continue
##'
##' @param Nmcmc positive integer; number of additional iterations to perform
##'
setMethod(
  "continue",
  signature=signature(object="pmcmcd_pomp"),
  function (object, Nmcmc = 1, ...) {

    ndone <- object@Nmcmc
    accepts <- object@accepts

    obj <- pmcmc(
      object=object,
      Nmcmc=Nmcmc,
      ...,
      .ndone=ndone,
      .accepts=accepts,
      .prev.pfp=as(object,"pfilterd_pomp"),
      .prev.log.prior=object@log.prior
    )

    obj@traces <- rbind(
      object@traces[,colnames(obj@traces)],
      obj@traces[-1,]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    ft <- array(dim=replace(dim(obj@filter.traj),2L,ndone+Nmcmc),
      dimnames=replace(dimnames(obj@filter.traj),2L,
        list(seq_len(ndone+Nmcmc))))
    ft[,seq_len(ndone),] <- object@filter.traj
    ft[,ndone+seq_len(Nmcmc),] <- obj@filter.traj
    obj@filter.traj <- ft
    obj@Nmcmc <- as.integer(ndone+Nmcmc)
    obj@accepts <- as.integer(accepts+obj@accepts)

    obj
  }
)
