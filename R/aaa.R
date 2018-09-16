##' @include pomp-package.R

.onAttach <- function (...) {
  packageStartupMessage("Welcome to pomp version 2\n",
    "For information on upgrading your pomp version < 2 code,\n",
    "see the 'pomp version 2 upgrade guide' at ",
    "https://kingaa.github.io/pomp/")
}
