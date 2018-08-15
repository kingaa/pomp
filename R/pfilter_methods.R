## methods for the 'pfilterd_pomp' class

setMethod(
  "cond.logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@cond.loglik
)

## 'coerce' method: allows for coercion of a 'pomp' object to a data-frame
setAs(
  from="pfilterd_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred.mean(from)
    pv <- pred.var(from)
    fm <- filter.mean(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      ess=eff.sample.size(from),
      cond.loglik=cond.logLik(from)
    )
    if (length(pm)>0)
      out <- cbind(out,pred.mean=t(pm))
    if (length(pv)>0)
      out <- cbind(out,pred.var=t(pv))
    if (length(fm)>0)
      out <- cbind(out,filter.mean=t(fm))
    out
  }
)

as.data.frame.pfilterd_pomp <- function (x, row.names, optional, ...) as(x,"data.frame")
