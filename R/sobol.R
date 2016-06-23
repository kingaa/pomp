sobol <- function (vars, n) {
    ep <- paste0("in ",sQuote("sobolDesign"),": ")
    if (!is.list(vars) || is.null(names(vars)))
        stop(ep,sQuote("vars")," must be a named list",call.=FALSE)
    if (!all(sapply(vars,function(x)is.numeric(x)&&(length(x)==2))))
        stop(ep,"each entry in ",sQuote("vars")," must specify a range",call.=FALSE)
    d <- length(vars)
    x <- tryCatch(
        .Call("sobol_sequence",as.integer(c(d,n))),
        error = function (e) {
            stop(ep,conditionMessage(e),call.=FALSE)
        }
    )
    y <- vapply(
        seq_len(d),
        function (k) {
            vars[[k]][1L]+(vars[[k]][2L]-vars[[k]][1L])*x[k,]
        },
        numeric(n)
    )
    colnames(y) <- names(vars)
    as.data.frame(y)
}

sobolDesign <- function (lower = numeric(0), upper = numeric(0), nseq) {
    ep <- paste0("in ",sQuote("sobolDesign"),": ")
    if (length(lower)!=length(upper))
        stop(ep,sQuote("lower")," and ",sQuote("upper")," must have same length",call.=FALSE)
    lnames <- names(lower)
    if (is.null(lnames))
        stop(ep,sQuote("lower")," and ",sQuote("upper")," must be named vectors",call.=FALSE)
    if (!all(sort(lnames)==sort(names(upper))))
        stop(ep,"names of ",sQuote("lower")," and ",sQuote("upper")," must match",call.=FALSE)
    upper <- upper[lnames]
    ranges <- lapply(lnames,function(k)c(lower[k],upper[k]))
    names(ranges) <- lnames
    sobol(ranges,n=as.integer(nseq))
}
