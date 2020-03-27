##' @include pomp-package.R
##' @importFrom utils packageVersion

.onAttach <- function (...) {
  packageStartupMessage(
    "Welcome to pomp version ",utils::packageVersion("pomp"),"!\n",
    "As of version 2.7.1.0, important changes have been made to the\n",
    "default settings of the particle filtering algorithms in\n",
    paste(sapply(c("pfilter","mif2","pmcmc","bsmc2"),sQuote),collapse=", "),".\n",
    "These changes are not backward compatible.\n",
    "See the package NEWS for the details.\n\n",
    "For information on upgrading your pomp version < 2 code, see the\n",
    dQuote("pomp version 2 upgrade guide")," at ",
    "https://kingaa.github.io/pomp/",".\n"
  )
}
