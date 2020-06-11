##' @include pomp-package.R

.onAttach <- function (...) {
  msg <- "
Welcome to pomp!
Version 3 incorporates some changes to the behavior of
package algorithms that are not backward compatible.
See the package NEWS for the details."
  packageStartupMessage("\n",paste(strwrap(msg),collapse="\n"),"\n")
}
