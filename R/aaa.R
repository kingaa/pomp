##' @include package.R

.onAttach <- function (...) {
  msg <- r"{
Welcome to pomp!

As of version 4.6, no user-visible pomp function has a name that
includes a dot (.). Function names have been changed to replace the dot
with an underscore (_). For more information, see the pomp blog:
https://kingaa.github.io/pomp/blog.html.}"
  packageStartupMessage("\n",paste(strwrap(msg),collapse="\n"),"\n")
}
