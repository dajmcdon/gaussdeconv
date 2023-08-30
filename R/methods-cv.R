#' @method summary cv_smoothed_mult
#' @export
summary.cv_smoothed_mult <- function(object, ...) {

  rlang::check_dots_empty()

  tab <- with(object, data.frame(
    lambda = lambda,
    index = seq_along(lambda),
    cv_scores = cv_scores,
    cv_se = cv_se,
    dof = full_fit$dof
  ))
  n <- nrow(tab)
  if (n > 5) {
    l1 <- which(abs(object$lambda - object$lambda.min) < 1e-10)
    l2 <- which(abs(object$lambda - object$lambda.1se) < 1e-10)
    idx <- c(1, l1, l2, n)
    tab <- tab[idx, ]
    rownames(tab) <- c("Min Lambda", "CV Minimizer", "1se Lambda", "Max Lambda")
  }

  out <- structure(
    list(call = object$call, table = tab, korder = object$full_fit$korder),
    class = "summary.cv_smoothed_mult")
  out
}

#' @method print cv_smoothed_mult
#' @export
print.cv_smoothed_mult <- function(x, ...) {
  print(summary(x, ...))
}


#' @method print summary.cv_smoothed_mult
#' @export
print.summary.cv_smoothed_mult <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...) {

  rlang::check_dots_empty()

  lambda_warning = NULL
  if (x$table$index[2] == 1) lambda_warning = "smallest"
  if (x$table$index[2] == x$table$index[4]) lambda_warning = "largest"

  cat("\nCall:", deparse(x$call), fill = TRUE)
  cat("\nDegree of the estimated piecewise polynomial curve:", x$korder, "\n")
  if (!is.null(lambda_warning)) {
    cat("Warning: the CV minimum occurred at the", lambda_warning,
        "lambda in the path.\n\n")
  }
  cat("\nSummary of cross validation across lambda:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' Plot cv_smoothed_mult
#'
#' @param x result of cv_estimate_smooth of class `cv_smoothed_mult`
#' @param which_lambda select which ct's to plot.
#'
#' If not provided, the cross validation score will be plotted. If provided a
#' list of lambda, the corresponding ct estimation will be plotted.
#'
#' If provided a string, it
#' must be either one of `lambda.min`, `lambda.1se`, or `cv_scores`.
#'
#'  * If provided `lambda.min`, plot ct which is generated from the lambda that
#'  minimizes the cross validation score.
#'
#'  * If provided `lambda.1se`, plot ct which is generated from the lambda whose
#'  corresponding cross validation score is 1 standard error away of the
#'  minimal cross validation score.
#'
#'  * If provided `cv_scores`, plot the cross validation score.
#'
#'  * If NULL, all estimated ct values are plotted.
#'
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_smooth(y, korder = 1, nfold = 3, nsol = 30)
#' plot(cv)
#' plot(cv, which_lambda = cv$lambda[1])
#' plot(cv, which_lambda = "lambda.min")
#' plot(cv, which_lambda = "lambda.1se")
#' plot(cv, NULL)
plot.cv_smoothed_mult <- function(
    x, which_lambda = c("cv_scores", "lambda.min", "lambda.1se"), ...) {

  rlang::check_dots_empty()
  plt_scores <- FALSE
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    if (which_lambda == "cv_scores") plt_scores <- TRUE
    else which_lambda <- x[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }

  if (plt_scores) {
    df <- with(x, data.frame(
      cv_scores = cv_scores,
      lambda = lambda,
      cv_se = cv_se,
      upper = cv_scores + cv_se,
      lower = cv_scores - cv_se
    ))
    plt <- ggplot2::ggplot(df) +
      ggplot2::geom_errorbar(ggplot2::aes(x = .data$lambda,
                                          y = .data$cv_scores,
                                          ymin = .data$lower,
                                          ymax = .data$upper,
                                          width = 0.1)) +
      ggplot2::geom_point(ggplot2::aes(x = .data$lambda, y = .data$cv_scores),
                          color = "darkblue") +
      ggplot2::geom_vline(xintercept = x$lambda.min, linetype = 'dotted') +
      ggplot2::geom_vline(xintercept = x$lambda.1se, linetype = 'dotted') +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Lambda", y = "CV scores") +
      ggplot2::scale_x_log10()
  } else {
    plt <- plot(x$full_fit, which_lambda)
  }

  return(plt)
}


#' Fitted cv_smoothed_mult
#'
#' @param object result of cross validation of type `cv_smoothed_mult`
#' @param which_lambda select which ct's to output. If not provided, all ct's
#' are returned. If provided a list of lambda, the corresponding ct estimate(s)
#' will be returned.
#'
#' If provided a string, it must be either one of `lambda.min` or `lambda.1se`.
#'
#'  * If provided `lambda.min`, return ct which is generated from
#'  the lambda that minimizes the cross validation score.
#'
#'  * If provided `lambda.1se`, return ct which is generated from the lambda
#'  whose corresponding cross validation score is 1 standard error away of the
#'  minimal cross validation score.
#' @param ... not used.
#'
#' @return ct's estimated from provided lambda
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_ct(y, korder = 3, nfold = 3, nsol = 30)
#' f <- fitted(cv)
#' f <- fitted(cv, which_lambda = cv$lambda[1])
#' f <- fitted(cv, which_lambda = "lambda.1se")
#' f <- fitted(cv, which_lambda = NULL)
fitted.cv_smoothed_mult <- function(
    object,
    which_lambda = c("lambda.min", "lambda.1se"),
    ...) {
  rlang::check_dots_empty()
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    which_lambda <- object[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }
  fitted(object$full_fit, which_lambda)
}


#' @importFrom stats coef
#' @export
coef.cv_smoothed_mult <- fitted.cv_smoothed_mult


#' Predict observed data using estimated Rt
#'
#' Given an object of class `cv_smoothed_mult` produced with [cv_estimate_ct()],
#' calculate predicted observations for the estimated ct values.
#' Note: This function is not intended for "new x" or to produce forecasts, but
#' rather to examine how ct relates to observables.
#'
#'
#' @param object result of cross validation of type `cv_smoothed_mult`
#' @param which_lambda Select which lambdas from the object to use. If not
#'   provided, all ct's are returned. Note that new lambdas not originally
#'   used in the estimation procedure may be provided, but the results will be
#'   calculated by linearly interpolating the estimated ct's.
#'
#'   The strings `lambda.min` or `lambda.1se` are allowed to choose either
#'   the lambda that minimizes the cross validation score or the largest lambda
#'   whose corresponding cross validation score is within 1 standard error of
#'   the minimal cross validation score.
#' @param ... not used.
#'
#' @return A vector or matrix of predicted case counts.
#' @export
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_ct(y, korder = 3, nfold = 3, nsol = 30)
#' p <- predict(cv)
#' p <- predict(cv, which_lambda = cv$lambda[1])
#' p <- predict(cv, which_lambda = "lambda.1se")
#' p <- predict(cv, which_lambda = NULL)
#' plot(y)
#' matlines(p, lty = 2)
predict.cv_smoothed_mult <- function(
    object,
    which_lambda = c("lambda.min", "lambda.1se"),
    ...) {
  rlang::check_dots_empty()
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    which_lambda <- object[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }
  predict(object$full_fit, which_lambda)
}

#' @export
interpolate_ct.cv_smoothed_mult <- function(
    object,
    xout,
    which_lambda = c("lambda.min", "lambda.1se"),
    ...) {
  rlang::check_dots_empty()
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    which_lambda <- object[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }
  interpolate_ct(object$full_fit, xout, lambda = which_lambda)
}
