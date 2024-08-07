##' C snippets
##'
##' Accelerating computations through inline snippets of C code
##'
##' \pkg{pomp} provides a facility whereby users can define their model's components using inline C code.
##' C snippets are written to a C file, by default located in the \R session's temporary directory, which is then compiled (via \code{\link[=SHLIB]{R CMD SHLIB}}) into a dynamically loadable shared object file.
##' This is then loaded as needed.
##'
##' @name Csnippet
##' @rdname csnippet
##' @include package.R
##' @family implementation information
##' @seealso spy
##'
##' @section Note to Windows and Mac users:
##'
##' By default, your \R installation may not support \code{\link[=SHLIB]{R CMD SHLIB}}.
##' The \href{https://kingaa.github.io/pomp/install.html}{package website contains installation instructions} that explain how to enable this powerful feature of \R.
##'
##' @inheritSection pomp Note for Windows users
##'
##' @section General rules for writing C snippets:
##'
##' In writing a C snippet one must bear in mind both the \emph{goal} of the snippet, i.e., what computation it is intended to perform, and the \emph{context} in which it will be executed.
##' These are explained here in the form of general rules.
##' Additional specific rules apply according to the function of the particular C snippet.
##' Illustrative examples are given in the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.
##' \enumerate{
##' \item C snippets must be valid C.
##' They will embedded verbatim in a template file which will then be compiled by a call to \code{\link[=SHLIB]{R CMD SHLIB}}.
##' If the resulting file does not compile, an error message will be generated.
##' Compiler messages will be displayed, but no attempt will be made by \pkg{pomp} to interpret them.
##' Typically, compilation errors are due to either invalid C syntax or undeclared variables.
##' \item State variables, parameters, observables, and covariates must be left undeclared within the snippet.
##' State variables and parameters are declared via the \code{statenames} or \code{paramnames} arguments to \code{pomp}, respectively.
##' Compiler errors that complain about undeclared state variables or parameters are usually due to failure to declare these in \code{statenames} or \code{paramnames}, as appropriate.
##' \item A C snippet can declare local variables.
##' Be careful not to use names that match those of state variables, observables, or parameters.
##' One must never declare state variables, observables, covariates, or parameters within a C snippet.
##' \item Names of observables must match the names given given in the data.
##' They must be referred to in measurement model C snippets (\code{rmeasure} and \code{dmeasure}) by those names.
##' \item If the \sQuote{pomp} object contains a table of covariates (see above), then the variables in the covariate table will be available, by their names, in the context within which the C snippet is executed.
##' \item Because the dot \sQuote{.} has syntactic meaning in C, \R variables with names containing dots (\sQuote{.}) are replaced in the C codes by variable names in which all dots have been replaced by underscores (\sQuote{_}).
##' \item The headers \file{R.h} and \file{Rmath.h}, provided with \R, will be included in the generated C file, making all of the \href{https://CRAN.R-project.org/doc/manuals/r-release/R-exts.html#The-R-API}{\R C API} available for use in the C snippet.
##' This makes a great many useful functions available, including all of \R's \href{https://CRAN.R-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions}{statistical distribution functions}.
##' \item The header \href{https://github.com/kingaa/pomp/blob/master/inst/include/pomp.h}{\file{pomp.h}}, provided with \pkg{pomp}, will also be included, making all of the \href{https://kingaa.github.io/pomp/C_API.html}{\pkg{pomp} C API} available for use in every C snippet.
##' \item Snippets of C code passed to the \code{globals} argument of \code{pomp} will be included at the head of the generated C file.
##' This can be used to declare global variables, define useful functions, and include arbitrary header files.
##' }
##'
##' @section Linking to precompiled libraries:
##' It is straightforward to link C snippets with precompiled C libraries.
##' To do so, one must make sure the library's header files are included;
##' the \code{globals} argument can be used for this purpose.
##' The \code{shlib.args} argument can then be used to specify additional arguments to be passed to \code{\link[=SHLIB]{R CMD SHLIB}}.
##' \href{https://kingaa.github.io/pomp/FAQ.html#linking-C-libraries}{FAQ 3.7} gives an example.
##'
##' @section C snippets are salted:
##' To prevent collisions in parallel computations, a \sQuote{pomp} object built using C snippets is \dQuote{salted} with the current time and a random number.
##' A result is that two \sQuote{pomp} objects, built on identical codes and data, will \strong{not} be identical as \R objects, though they will be functionally identical in every respect.
##'
NULL

setClass(
  "Csnippet",
  slots=c(
    text="character"
  ),
  prototype=prototype(
    text=character(0)
  )
)

##' @name Csnippet
##' @rdname csnippet
##' @param text character; text written in the C language
##' @export
Csnippet <- function (text) {
  new("Csnippet",text=as.character(text))
}

setAs(
  from="Csnippet",
  to="character",
  def = function (from) {
    from@text
  }
)

##' @export
as.character.Csnippet <- function(x, ...) as(x,"character")
