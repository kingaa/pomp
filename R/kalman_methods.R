## methods for the 'kalmand_pomp' class

setMethod(
  "cond.logLik",
  signature=signature(object="kalmand_pomp"),
  definition=function(object,...)object@cond.loglik
)

setAs(
  from="kalmand_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred.mean(from)
    fm <- filter.mean(from)
    fc <- forecast(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      cond.loglik=cond.logLik(from)
    )
    if (length(pm)>0)
      out <- cbind(out,pred.mean=t(pm))
    if (length(fm)>0)
      out <- cbind(out,filter.mean=t(fm))
    if (length(fc)>0)
      out <- cbind(out,forecast=t(fc))
    out
  }
)

as.data.frame.kalmand_pomp <- function (x, row.names, optional, ...) {
  as(x,"data.frame")
}
