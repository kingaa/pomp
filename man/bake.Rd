\name{bake}
\alias{bake}
\alias{stew}
\title{Run a computation only if needed, saving the results for later use.}
\description{
  On cooking shows, recipes requiring lengthy baking or stewing are prepared beforehand.
  The \code{bake} and \code{stew} functions perform analogously:
  an \R computation is performed and stored in a named file.
  If the function is called again and the file is present, the computation is not executed; rather, the results are loaded from the file in which they were previously stored.
}
\usage{
bake(file, expr)
stew(file, expr)
}
\arguments{
  \item{file}{
    Name of the binary data file in which the result will be stored or retrieved, as appropriate.
    For \code{bake}, this will contain a single \R object and hence be an RDS file (extension \sQuote{rds});
    for \code{stew}, this will contain one or more named \R objects and hence be an RDA file (extension \sQuote{rda}).
  }
  \item{expr}{
    Expression to be evaluated.
  }
}
\examples{
  \dontrun{
  bake(file="example1.rds",{
    x <- runif(1000)
    mean(x)
  })

  stew(file="example2.rds",{
    x <- runif(10)
    y <- rnorm(n=10,mean=3*x+5,sd=2)
  })

  plot(x,y)  
  }
}
\details{
  Both \code{bake} and \code{stew} first test to see whether \code{file} exists.
  If it does, \code{bake} reads it using \code{\link{readRDS}} and returns the resulting object.
  By contrast, \code{stew} loads the file using \code{\link{load}} and copies the objects it contains into the user's workspace (or the environment of the call to \code{stew}).

  If \code{file} does not exist, then both \code{bake} and \code{stew} evaluate the expression \code{expr};
  they differ in the results that they save.
  \code{bake} saves the value of the evaluated expression to \code{file} as a single \R object.
  The name of that object is not saved.
  By contrast, \code{stew} creates a local environment within which \code{expr}is evaluated;
  all objects in that environment are saved (by name) in \code{file}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
