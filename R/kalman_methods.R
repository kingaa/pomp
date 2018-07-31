## methods for the 'kalmand.pomp' class

setMethod(
    "logLik",
    signature=signature(object="kalmand.pomp"),
    definition=function(object,...)object@loglik
)

setMethod(
    "cond.logLik",
    signature=signature(object="kalmand.pomp"),
    definition=function(object,...)object@cond.loglik
)

## 'coerce' method: allows for coercion of a 'pomp' object to a data-frame
setAs(
    from="kalmand.pomp",
    to="data.frame",
    def = function (from) {
        pm <- pred.mean(from)
        fm <- filter.mean(from)
        out <- cbind(
            as(as(from,"pomp"),"data.frame"),
            cond.loglik=cond.logLik(from)
        )
        if (length(pm)>0)
            out <- cbind(out,pred.mean=t(pm))
        if (length(fm)>0)
            out <- cbind(out,filter.mean=t(fm))
        out
    }
)

as.data.frame.kalmand.pomp <- function (x, row.names, optional, ...) {
    as(x,"data.frame")
}

## extract the prediction means
setMethod(
    "pred.mean",
    signature=signature(object="kalmand.pomp"),
    definition=function (object, pars, ...) {
        if (missing(pars)) pars <- rownames(object@pred.mean)
        object@pred.mean[pars,,drop=FALSE]
    }
)

## extract the filtering means
setMethod(
    "filter.mean",
    signature=signature(object="kalmand.pomp"),
    definition=function (object, pars, ...) {
        if (missing(pars)) pars <- rownames(object@filter.mean)
        object@filter.mean[pars,,drop=FALSE]
    }
)
