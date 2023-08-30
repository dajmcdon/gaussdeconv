#' @method summary smoothed_mult
#' @export
summary.smoothed_mult <- function(object, ...) {
  rlang::check_dots_empty()
  ns <- length(object$lambda)
  if (ns > 5) {
    xlam <- round(stats::quantile(1:ns))
    names(xlam) <- rev(c("Max.", "3rd Qu.", "Median", "1st Qu.", "Min."))
  } else {
    xlam <- seq_len(ns)
    names(xlam) <- paste0("s", seq_len(ns))
  }
  tab <- with(object, data.frame(
    lambda = lambda[xlam],
    index = xlam,
    approx_dof = dof[xlam],
    niterations = niter[xlam]))
  lambda <- object$lambda
  rownames(tab) <- names(xlam)
  out <- structure(
    list(call = object$call, table = tab, korder = object$korder, nlam = ns),
    class = "summary.smoothed_mult")
  out
}

#' @method print summary.smoothed_mult
#' @export
print.summary.smoothed_mult <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...) {
  rlang::check_dots_empty()
  cat("\nCall: ", deparse(x$call), fill = TRUE)
  cat("\nDegree of the estimated piecewise polynomial curve:", x$korder, "\n")
  cat("\nSummary of the", x$nlam, "estimated models:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' @method print smoothed_mult
#' @export
print.smoothed_mult <- function(
    x,
    digits = min(3, getOption("digits") - 3),
    ...) {
  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}

#' Plot estimated ct values from a `smoothed_mult` object
#'
#' Produces a figure showing some or all estimated Rt values for different
#' values of the penalty. The result is a [ggplot2::ggplot()]. Additional user
#' modifications can be added as desired.
#'
#'
#' @param x output of the function [estimate_ct()] of class `smoothed_mult`
#' @param lambda select which ct's to plot. If not provided,
#'   all ct's are plotted.
#' @param ... Not used.
#'
#' @export
#'
#' @importFrom rlang .data
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_ct(y, lambda = log(c(1.1, 1.3, 1.5)))
#' plot(out)
plot.smoothed_mult <- function(x, lambda = NULL, ...) {
  arg_is_positive(lambda, allow_null = TRUE)

  n <- length(x$y)
  if (is.null(lambda)) {
    ct <- x$ct
    lambda <- x$lambda
  } else {
    ct <- fitted(x, lambda = lambda)
  }
  k <- length(lambda)

  df <- data.frame(
    ct = c(ct),
    lambda = rep(lambda, each = n),
    Time = rep(x$x, k)
  )

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$Time,
      y = .data$ct,
      colour = .data$lambda,
      group = .data$lambda)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_viridis_c(trans = "log10")
  if (k == 1) plt <- plt + ggplot2::theme(legend.position = "none")
  plt
}

#' @importFrom stats fitted
#' @export
fitted.smoothed_mult <- function(object, lambda = NULL, ...) {
  rlang::check_dots_empty()
  arg_is_positive(lambda, allow_null = TRUE)

  if (is.null(lambda)) return(object$ct)

  lam_list <- interpolate_lambda(object$lambda, lambda)
  k <- length(lambda)
  ct <- object$ct
  ret <- ct[ ,lam_list$left, drop = FALSE] %*% diag(lam_list$frac, k, k) +
    ct[ ,lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, k, k)
  drop(ret)
}

#' @importFrom stats coef
#' @export
coef.smoothed_mult <- fitted.smoothed_mult

#' Predict observed data using estimated ct
#'
#' Given an object of class `smoothed_mult` produced with [estimate_ct()],
#' calculate predicted observed cases for the estimated ct values.
#' Note: This function is not intended for "new x" or to produce forecasts, but
#' rather to examine how ct relates to observables.
#'
#' @param object An object of class `smoothed_mult` produced with [estimate_ct()].
#' @param lambda Select which lambdas from the object to use. If not provided
#'   (the default), all are returned. Note that new lambdas not originally
#'   used in the estimation procedure may be provided, but the results will be
#'   calculated by linearly interpolating the estimated ct's.
#' @param ... Not used.
#'
#' @return A vector or matrix of predicted case counts.
#' @export
#'
#' @importFrom stats predict
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_ct(y, nsol = 10)
#' preds <- predict(out)
#' plot(y)
#' matlines(preds, lty = 1)
predict.smoothed_mult <- function(object, lambda = NULL, ...) {
  rlang::check_dots_empty()

  ct <- fitted(object, lambda)
  ct * object$convolved
}

#' @export
interpolate_ct.smoothed_mult <- function(object, xout, lambda = NULL, ...) {
  rlang::check_dots_empty()
  xin <- object$x
  if (inherits(xin, "Date")) xin <- as.numeric(xin)
  arg_is_positive(lambda, allow_null = TRUE)
  if (is.unsorted(xout)) xout <- sort(xout)

  ct <- fitted(object, lambda = lambda)
  interp <- apply(ct, 2, function(zz) {
    dspline::dspline_interp(zz, object$korder, xin, xout)
  })
  interp
}

