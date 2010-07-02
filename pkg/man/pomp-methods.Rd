\name{pomp-methods}
\docType{methods}
\alias{pomp-methods}
\alias{plot,pomp-method}
\alias{plot-pomp}
\alias{print,pomp-method}
\alias{print-pomp}
\alias{show,pomp-method}
\alias{show-pomp}
\alias{time,pomp-method}
\alias{time-pomp}
\alias{time<-}
\alias{time<-,pomp-method}
\alias{time<--pomp}
\alias{data.array,pomp-method}
\alias{data.array-pomp}
\alias{data.array}
\alias{states,pomp-method}
\alias{states-pomp}
\alias{states}
\alias{as,pomp-method}
\alias{coerce,pomp,data.frame-method}
\alias{coef,pomp-method}
\alias{coef-pomp}
\alias{coef<-}
\alias{coef<-,pomp-method}
\alias{coef<--pomp}
\title{Methods of the "pomp" class}
\description{Methods of the \code{pomp} class.}
\usage{
\S4method{coef}{pomp}(object, pars, \dots)
\S4method{coef}{pomp}(object, pars, \dots) <- value
\S4method{data.array}{pomp}(object, vars, \dots)
\S4method{states}{pomp}(object, vars, \dots)
\S4method{time}{pomp}(x, t0 = FALSE, \dots)
\S4method{time}{pomp}(object, t0 = FALSE, \dots) <- value
\S4method{show}{pomp}(object)
\S4method{as}{pomp}(object, class)
\S4method{coerce}{pomp,data.frame}(from, to = "data.frame", strict = TRUE)
\S4method{print}{pomp}(x, \dots)
\S4method{plot}{pomp}(x, y, variables, panel = lines,
      nc = NULL, yax.flip = FALSE,
      mar = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
      oma = c(6, 0, 5, 0), axes = TRUE, \dots)
}
\arguments{
  \item{object, x}{The \code{pomp} object.}
  \item{pars}{
    optional character;
    names of parameters to be retrieved or set.
  }
  \item{vars}{
    optional character;
    names of observed variables to be retrieved.
  }
  \item{value}{
    numeric;
    values to be assigned to the parameters.
  }
  \item{t0}{
    logical;
    if TRUE on a call to \code{time}, the zero time is prepended to the time vector;
    if TRUE on a call to \code{time<-}, the first element in \code{value} is taken to be the initial time.
  }
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{from, to}{
    the classes between which coercion should be performed.
  }
  \item{strict}{
    ignored.
  }
  \item{y}{ignored.}
  \item{variables}{
    optional character;
    names of variables to plot.
  }
  \item{panel}{
    a function of prototype \code{panel(x, col, bg, pch, type, ...)} which gives the action to be carried out in each panel of the display.
  }
  \item{nc}{
    the number of columns to use.
    Defaults to 1 for up to 4 series, otherwise to 2.
  }
  \item{yax.flip}{
    logical;
    if TRUE, the y-axis (ticks and numbering) should flip from side 2 (left) to 4 (right) from series to series.
  }
  \item{mar, oma}{
    the 'par' settings for 'mar' and 'oma' to use.
    Modify with care!
  }
  \item{axes}{
    logical;
    indicates if x- and y- axes should be drawn.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\details{
  \describe{
    \item{coef}{
      \code{coef(object)} returns the contents of the \code{params} slot of \code{object}.
      \code{coef(object,pars)} returns only those parameters named in \code{pars}.
    }
    \item{coef<-}{
      Assigns values to the \code{params} slot of the \code{pomp} object.
      If \code{coef(object)} exists, then \code{coef(object) <- value} has the effect of replacing the parameters of \code{object} with \code{value};
      the names of \code{value} will be ignored and the names of \code{coef(object)} will be unchanged.
      If \code{coef(object)} does not exist, then \code{coef(object) <- value} assigns \code{value} to the parameters of \code{object};
      the names of \code{coef(object)} will be those of \code{value} and an error will be generated if \code{value} does not have names.
      If \code{coef(object)} exists, then \code{coef(object,pars) <- value} replaces those parameters of \code{object} named in \code{pars} with the elements of \code{value};
      the names of \code{value} are ignored.
      If \code{coef(object)} does not exist, then \code{coef(object,pars) <- value} assigns \code{value} to the parameters of \code{object};
      in this case, the names of \code{object} will be \code{pars} and the names of \code{value} will be ignored.
    }
    \item{data.array}{
      \code{data.array(object)} returns the array of observations.
      \code{data.array(object,vars)} gives just the observations of variables named \code{vars}.
      \code{vars} may specify the variables by position or by name.
    }
    \item{states}{
      \code{states(object)} returns the array of states.
      \code{states(object,vars)} gives just the state variables named in \code{vars}.
      \code{vars} may specify the variables by position or by name.
    }
    \item{time}{
      \code{time(object)} returns the vector of observation times.
      \code{time(object,t0=TRUE)} returns the vector of observation times with the zero-time \code{t0} prepended.
    }
    \item{time<-}{
      \code{time(object) <- value} replaces the observation times slot (\code{times}) of \code{object} with \code{value}.
      \code{time(object,t0=TRUE) <- value} has the same effect, but the first element in \code{value} is taken to be the initial time.
      The second and subsequent elements of \code{value} are taken to be the observation times.
      Those data and states (if they exist) corresponding to the new times are retained.
    }
    \item{show}{Displays the \code{pomp} object.}
    \item{plot}{
      Plots the data and state trajectories (if the latter exist).
      Additional arguments are passed to the low-level plotting routine.
    }
    \item{print}{Prints the \code{pomp} object in a nice way.}
    \item{as, coerce}{
      The \code{coerce} method should typically not be used directly.
      It is defined by \code{setAs} as a method to be used by \code{as}.
      A \code{pomp} object can be coerced to a data frame via \code{as(object,"data.frame")}.
      The data frame contains the times, the data, and the state trajectories, if they exist.
    }
    \item{rprocess}{
      simulates the process model.
      See \code{\link{rprocess-pomp}}.
    }
    \item{dprocess}{
      evaluates the process model density.
      See \code{\link{dprocess-pomp}}.
    }
    \item{rmeasure}{
      simulates the measurement model.
      See \code{\link{rmeasure-pomp}}.
    }
    \item{dmeasure}{
      evaluates the measurement-model density.
      See \code{\link{dmeasure-pomp}}.
    }
    \item{skeleton}{
      evaluates the deterministic skeleton (be it a vector field or a map).
      See \code{\link{skeleton-pomp}}.
    }
    \item{init.state}{
      returns a vector of initialial conditions.
      See \code{\link{init.state-pomp}}.
    }
    \item{simulate}{
      \code{simulate} can be used to simulate state and observation trajectories.
      See documentation under \code{\link{simulate-pomp}}.
    }
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}},
  \link{pomp-class},
  \code{\link{rprocess}},
  \code{\link{dprocess}},
  \code{\link{rmeasure}},
  \code{\link{dmeasure}},
  \code{\link{init.state}},
  \link{simulate-pomp}
}
\keyword{models}
\keyword{ts}
