#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @inheritParams estimate_ct
#' @param nfold Integer. This number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds. Smallest allowable value is `nfold = 2`.
#' @param error_measure Metric used to calculate cross validation scores.
#'   Must be choosen from `mse` or `mae`.
#'   `mse` calculates the mean square error;
#'   `mae` calculates the mean absolute error;
#' @param ... additional parameters passed to `estimate_ct()` function

#' @return An object with S3 class `"cv_smoothed_mult"`. Among the list components:
#' * `full_fit` An object with S3 class `"smoothed_mult"`, fitted with all
#' `observed_counts` and `lambda`
#' * `cv_scores` leave-kth-out cross validation scores
#' * `cv_se` leave-kth-out cross validation standard error
#' * `lambda.min` lambda which achieved the optimal cross validation score
#' * `lambda.1se` lambda that gives the optimal cross validation score
#' within one standard error.
#' * `lambda` the value of `lambda` used in the algorithm.
#' @export
#'
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_ct(y, korder = 3, nfold = 3, nsol = 30)
#' cv
cv_estimate_ct <- function(
    y,
    korder = 3L,
    dist_gamma = c(2.5, 2.5),
    nfold = 3L,
    error_measure = c("mse", "mae"),
    x = 1:n,
    lambda = NULL,
    maxiter = 1e6L,
    delay_distn = NULL,
    ...) {

  arg_is_pos_int(nfold)
  n <- length(y)
  arg_is_length(n, x)
  xin <- x
  if (inherits(xin, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)

  if (nfold == 1) cli_abort("nfold must be greater than 1")

  ## Run program one time to create lambda
  full_data_fit <- estimate_ct(
    y = y,
    korder = korder,
    dist_gamma = dist_gamma,
    x = xin,
    lambda = lambda,
    maxiter = maxiter,
    delay_distn = delay_distn,
    ...)

  if (is.null(lambda)) lambda <- full_data_fit$lambda

  foldid <- c(0, rep_len(1:nfold, n - 2), 0)
  cvall <- matrix(0, nfold, length(lambda))
  error_measure <- match.arg(error_measure)
  err_fun <- switch(
    error_measure,
    mse = function(y, m) (y - m)^2,
    mae = function(y, m) abs(y - m),
    deviance = function(y, m) (y - m)^2 / 2
  )

  for (f in 1:nfold) {
    train_idx <- foldid != f
    test_idx <- foldid == f

    mod <- estimate_ct(
      y = y[train_idx],
      korder = korder,
      dist_gamma = dist_gamma,
      x = x[train_idx],
      lambda = lambda,
      maxiter = maxiter,
      delay_distn = delay_distn,
      ...)

    interp_ct <- interpolate_ct(mod, x[test_idx])

    wpc <- delay_calculator(
      y = y[train_idx],
      x = x[train_idx] - min(x[train_idx]) + 1,
      dist_gamma = dist_gamma,
      delay_distn = delay_distn,
      output_partial_seq = FALSE)


    preds <- interp_ct * wpc[test_idx]
    score <- colMeans(err_fun(y[test_idx], preds))
    cvall[f,] <- score
  }

  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2) / sqrt(nfold)
  i0 <- which.min(cv_scores)

  structure(
    list(
      full_fit = full_data_fit,
      cv_scores = cv_scores,
      cv_se = cv_se,
      lambda = lambda,
      lambda.min = lambda[which.min(cv_scores)],
      lambda.1se = max(lambda[cv_scores <= cv_scores[i0] + cv_se[i0]]),
      call = match.call()
    ),
    class = "cv_smoothed_mult"
  )
}
