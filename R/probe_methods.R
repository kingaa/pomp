setAs(
  from="probed_pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(rbind(from@datvals,from@simvals))
    row.names(x) <- seq.int(from=0,to=nrow(x)-1)
    x$.id <- factor(c("data",rep("sim",nrow(x)-1)))
    x
  }
)

setMethod(
  "probevals",
  signature=signature(object="probed_pomp"),
  definition=function (object, ...) {
    dv <- object@datvals
    list(
      data=array(dv,dim=c(1,length(dv)),
        dimnames=list(rep="data",probe=names(dv))),
      sims=object@simvals
    )
  }
)


as.data.frame.probed_pomp <- function (x, row.names, optional, ...)
  as(x,"data.frame")
