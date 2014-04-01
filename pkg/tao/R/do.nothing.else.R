do.nothing.else <- function (...) {
  while (TRUE) {
    do.nothing(...)
  }
  invisible(NULL)
}
