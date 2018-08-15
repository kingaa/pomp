## methods to concatenate objects into useful listies
setGeneric(
    "concat",
    function (...)
        standardGeneric("concat")
)

setMethod(
  "concat",
  signature=signature(...="missing"),
  definition=function(...) {
      NULL
  }
)

setMethod(
  "concat",
  signature=signature(...="ANY"),
  definition=function(...) {
    stop(sQuote("c")," is not defined for objects of mixed class.",
      call.=FALSE)
  }
)

setMethod(
  "concat",
  signature=signature(...="Pomp"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pompList",unlist(y))
  }
)

setMethod(
  "concat",
  signature=signature(...="Pfilter"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pfilterList",unlist(y))
  }
)

setMethod(
  "concat",
  signature=signature(...="Abc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("abcList",unlist(y))
  }
)

setMethod(
  "concat",
  signature=signature(...="Mif2"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("mif2List",unlist(y))
  }
)

setMethod(
  "concat",
  signature=signature(...="Pmcmc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pmcmcList",unlist(y))
  }
)

c.Pomp <- concat
c.Pfilter <- concat
c.Abc <- concat
c.Mif2 <- concat
c.Pmcmc <- concat
