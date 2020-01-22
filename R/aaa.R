##' @include pomp-package.R

.onAttach <- function (...) {
  packageStartupMessage(
    "Welcome to pomp version 2!\n",
    "For information on upgrading your pomp version < 2 code, see the\n",
    sQuote("pomp version 2 upgrade guide")," at ",
    "https://kingaa.github.io/pomp/",".",
    "\n\nAlso, note that the default value of the ",sQuote("tol"),
    " argument of ",sQuote("pfilter"),"\n(and other particle-filter functions, i.e., ",
    sQuote("mif2"),", ",sQuote("pmcmc"),", ",sQuote("bsmc2"),") has changed to 0.\n",
    "The ",sQuote("tol")," argument is now deprecated and will be removed in a future release."
  )
}
