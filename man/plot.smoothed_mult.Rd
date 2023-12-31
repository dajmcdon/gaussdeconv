% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{plot.smoothed_mult}
\alias{plot.smoothed_mult}
\title{Plot estimated ct values from a \code{smoothed_mult} object}
\usage{
\method{plot}{smoothed_mult}(x, lambda = NULL, ...)
}
\arguments{
\item{x}{output of the function \code{\link[=estimate_ct]{estimate_ct()}} of class \code{smoothed_mult}}

\item{lambda}{select which ct's to plot. If not provided,
all ct's are plotted.}

\item{...}{Not used.}
}
\description{
Produces a figure showing some or all estimated Rt values for different
values of the penalty. The result is a \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}}. Additional user
modifications can be added as desired.
}
\examples{
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_ct(y, lambda = log(c(1.1, 1.3, 1.5)))
plot(out)
}
