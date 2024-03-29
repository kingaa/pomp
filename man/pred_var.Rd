% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_var.R
\name{pred_var}
\alias{pred_var}
\alias{pred_var,ANY-method}
\alias{pred_var,missing-method}
\alias{pred_var,pfilterd_pomp-method}
\title{Prediction variance}
\usage{
\S4method{pred_var}{pfilterd_pomp}(object, vars, ..., format = c("array", "data.frame"))
}
\arguments{
\item{object}{result of a filtering computation}

\item{vars}{optional character; names of variables}

\item{...}{ignored}

\item{format}{format of the returned object}
}
\description{
The variance of the prediction distribution
}
\details{
The prediction distribution is that of
\deqn{X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1},}{Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*,}
where \eqn{X(t_k)}{Xk}, \eqn{Y(t_k)}{Yk} are the latent state and observable processes, respectively, and \eqn{y^*_k}{yk*} is the data, at time \eqn{t_k}{tk}.

The prediction variance is therefore the variance of this distribution
\deqn{\mathrm{Var}[X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1}].}{Var[Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*].}
}
\seealso{
More on sequential Monte Carlo methods: 
\code{\link{bsmc2}()},
\code{\link{cond_logLik}()},
\code{\link{eff_sample_size}()},
\code{\link{filter_mean}()},
\code{\link{filter_traj}()},
\code{\link{kalman}},
\code{\link{mif2}()},
\code{\link{pfilter}()},
\code{\link{pmcmc}()},
\code{\link{pred_mean}()},
\code{\link{saved_states}()},
\code{\link{wpfilter}()}

Other extraction methods: 
\code{\link{coef}()},
\code{\link{cond_logLik}()},
\code{\link{covmat}()},
\code{\link{eff_sample_size}()},
\code{\link{filter_mean}()},
\code{\link{filter_traj}()},
\code{\link{forecast}()},
\code{\link{logLik}},
\code{\link{obs}()},
\code{\link{pred_mean}()},
\code{\link{saved_states}()},
\code{\link{spy}()},
\code{\link{states}()},
\code{\link{summary}()},
\code{\link{time}()},
\code{\link{timezero}()},
\code{\link{traces}()}
}
\concept{extraction methods}
\concept{particle filter methods}
