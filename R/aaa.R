##' @include pomp-package.R

.onAttach <- function (...) {
  packageStartupMessage(
    "Welcome to pomp!\n",
    "Version 3 incorporates some changes to the behavior of\n",
    "package algorithms that are not backward compatible.\n",
    "See the package NEWS for the details.\n\n"
  )
}
