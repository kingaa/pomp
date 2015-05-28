\name{pomp methods}
\docType{methods}
\alias{pomp methods}
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
\alias{timezero}
\alias{timezero,pomp-method}
\alias{timezero-pomp}
\alias{timezero<-}
\alias{timezero<-,pomp-method}
\alias{timezero<--pomp}
\alias{window,pomp-method}
\alias{window-pomp}
\alias{obs,pomp-method}
\alias{obs-pomp}
\alias{obs}
\alias{data.array,pomp-method}
\alias{data.array-pomp}
\alias{data.array}
\alias{states,pomp-method}
\alias{states-pomp}
\alias{states}
\alias{as,pomp-method}
\alias{coerce,pomp,data.frame-method}
\alias{as.data.frame.pomp}
\alias{partrans}
\alias{partrans,pomp-method}
\alias{partrans-pomp}
\alias{coef,pomp-method}
\alias{coef-pomp}
\alias{coef<-}
\alias{coef<-,pomp-method}
\alias{coef<--pomp}
\title{Functions for manipulating, displaying, and extracting information from objects of the \code{pomp} class}
\description{This page documents the various methods that allow one to extract information from, display, plot, and modify \code{pomp} objects.}
\usage{
\S4method{coef}{pomp}(object, pars, transform = FALSE, \dots)
\S4method{coef}{pomp}(object, pars, transform = FALSE, \dots) <- value
\S4method{obs}{pomp}(object, vars, \dots)
\S4method{partrans}{pomp}(object, params, dir = c("fromEstimationScale",
      "toEstimationScale","forward","inverse"), \dots)
\S4method{plot}{pomp}(x, y, variables, panel = lines,
      nc = NULL, yax.flip = FALSE,
      mar = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
      oma = c(6, 0, 5, 0), axes = TRUE, \dots)
\S4method{print}{pomp}(x, \dots)
\S4method{show}{pomp}(object)
\S4method{states}{pomp}(object, vars, \dots)
\S4method{time}{pomp}(x, t0 = FALSE, \dots)
\S4method{time}{pomp}(object, t0 = FALSE, \dots) <- value
\S4method{timezero}{pomp}(object, \dots)
\S4method{timezero}{pomp}(object, \dots) <- value
\S4method{window}{pomp}(x, start, end, \dots)
\S4method{as}{pomp}(object, class)
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
  \item{transform}{
    optional logical;
    should the parameter transformations be applied?
  }
  \item{value}{
    numeric;
    values to be assigned.
  }
  \item{params}{
    a vector or matrix of parameters to be transformed.
  }
  \item{dir}{
    direction of the transformation.
    \code{dir="forward"} applies the transformation from the \dQuote{natural} scale to the \dQuote{internal} scale.
    This is the transformation specified by the \code{parameter.transform} argument to \code{pomp};
    it is stored in the \sQuote{par.trans} slot of \code{object}.
    \code{dir="inverse"} applies the inverse transformation (stored in the \sQuote{par.untrans} slot).
  }
  \item{t0}{
    logical;
    if TRUE on a call to \code{time}, the zero time is prepended to the time vector;
    if TRUE on a call to \code{time<-}, the first element in \code{value} is taken to be the initial time.
  }
  \item{start, end}{
    start and end times of the window.
  }
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
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
      \preformatted{coef(object,transform=TRUE)} returns \preformatted{parameter.inv.transform(coef(object))}, where \code{parameter.inv.transform} is the user parameter inverse transformation function specified when \code{object} was created.
      Likewise, \preformatted{coef(object,pars,transform=TRUE)} returns \preformatted{parameter.inv.transform(coef(object))[pars]}.
    }
    \item{coef<-}{
      Assigns values to the \code{params} slot of the \code{pomp} object.
      \code{coef(object) <- value} has the effect of replacing the parameters of \code{object} with \code{value}.
      If \code{coef(object)} exists, then \code{coef(object,pars) <- value} replaces those parameters of \code{object} named in \code{pars} with the elements of \code{value};
      the names of \code{value} are ignored.
      If some of the names in \code{pars} do not already name parameters in \code{coef(object)}, then they are concatenated.
      If \code{coef(object)} does not exist, then \code{coef(object,pars) <- value} assigns \code{value} to the parameters of \code{object};
      in this case, the names of \code{object} will be \code{pars} and the names of \code{value} will be ignored.
      \code{coef(object,transform=TRUE) <- value} assigns \code{parameter.transform(value)} to the \code{params} slot of \code{object}.
      Here, \code{parameter.transform} is the parameter transformation function specified when \code{object} was created.
      \code{coef(object,pars,transform=TRUE) <- value} first, discards any names the \code{value} may have, sets \code{names(value) <- pars}, and then replaces the elements of \code{object}'s \code{params} slot \code{parameter.transform(value)}.
      In this case, if some of the names in \code{pars} do not already name parameters in \code{coef(object,transform=TRUE)}, then they are concatenated.
    }
    \item{obs}{
      \code{obs(object)} returns the array of observations.
      \code{obs(object,vars)} gives just the observations of variables named in \code{vars}.
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
    \item{timezero, timezero<-}{
      \code{timezero(object)} returns the zero-time \code{t0}.
      \code{timezero(object) <- value} sets the zero-time to \code{value}.
    }
    \item{window}{
      \code{window(x,start=t1,end=t2} returns a new \code{pomp} object, identical to \code{x} but with only the data in the window between times \code{t1} and \code{t2} (inclusive).
      By default, \code{start} is the time of the first observation and \code{end} is the time of the last.
    }
    \item{show}{Displays the \code{pomp} object.}
    \item{plot}{
      Plots the data and state trajectories (if the latter exist).
      Additional arguments are passed to the low-level plotting routine.
    }
    \item{print}{Prints the \code{pomp} object in a nice way.}
    \item{as}{
      A \code{pomp} object can be coerced to a data frame via \preformatted{as(object,"data.frame").}
      The data frame contains the times, the data, and the state trajectories, if they exist.
    }
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}},
  \link{pomp low-level interface},
  \link[=simulate-pomp]{simulate},
  \code{\link{pfilter}}, \code{\link{probe}}.
}
\keyword{programming}
\keyword{ts}
