##' Compartmental epidemiological models
##'
##' @description
##' \code{sir} is a \sQuote{pomp} object encoding a simple seasonal SIR model.
##' Simulation is performed using an Euler multinomial approximation.
##'
##' \code{sir2} has the same model implemented using Gillespie's algorithm.
##'
##' \code{bbs} is a nonseasonal SIR model together with data from a 1978
##' outbreak of influenza in a British boarding school.
##'
##' This and similar examples are discussed and constructed in tutorials
##' available on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' The boarding school influenza outbreak is described in Anonymous (1978).
##'
##' @details
##' Do \preformatted{
##'     pompExample(sir,show=TRUE)
##'     pompExample(sir2,show=TRUE)
##'     pompExample(bbs,show=TRUE)
##' }
##' to see the codes that generate these examples.
##'
##' @name sir_models
##' @rdname sir
##' @aliases sir sir2 bbs
##' @docType data
##' @references Anonymous (1978).  Influenza in a boarding school.  British
##' Medical Journal 1:587
##' @keywords datasets models
##' @family pomp examples
##' @examples
##'
##' pompExample(sir)
##' plot(sir)
##' plot(simulate(sir))
##' coef(sir)
##'
##' pompExample(sir2)
##' plot(sir2)
##' plot(simulate(window(sir2,end=3)))
##' coef(sir2)
##'
##' pompExample(bbs)
##' plot(bbs)
##' coef(bbs)
##' as.data.frame(bbs)
##'
NULL
