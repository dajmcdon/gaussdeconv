#' Estimate smoothed ratio between observations and weighted past
#'
#'
#' This function solves the smoothness penalized Gaussian
#' regression (trend filtering) of the form:
#' \deqn{\hat{\theta} = \arg\min_{\theta} \frac{1}{n} \sum_{i=1}^n (y_i -
#' \theta_i w_i)^2 + \lambda||D^{(k+1)}\theta||_1, }
#' where \eqn{y_i} are the observations at day
#' \eqn{i}, \eqn{w_i} are the weighted past observations at day \eqn{i},
#' convolved with a delay distribution (a gamma density by default),
#' \eqn{\lambda} is a tuning parameter, larger values resulting in smoother
#' estimates, and \eqn{D^{(k+1)}} is the \eqn{(k+1)}-th order
#' difference matrix.
#'
#' @param y vector of the observations
#' @param korder Integer. Degree of the piecewise polynomial curve to be
#'   estimated. For example, `korder = 0` corresponds to a piecewise constant
#'   curve.
#' @param lambda Vector. A user supplied sequence of tuning parameters which
#'   determines the balance between data fidelity and
#'   smoothness of the estimated ratio; larger `lambda` results in a smoother
#'   estimate. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the largest value
#'   of `lambda` that would result in a maximally smooth estimate, and `lambda_min_ratio`.
#'   Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param maxiter Integer. Maximum number of iterations for the estimation
#'   algorithm.
#' @param init a list of internal configuration parameters of class
#'   `rt_admm_configuration`.
#' @param dist_gamma Vector of length 2. These are the shape and scale for the
#'   assumed delay distribution.
#' @param x a vector of positions at which the observations have occurred. In an
#'   ideal case, we would observe data at regular intervals (e.g. daily or
#'   weekly) but this may not always be the case. May be numeric or Date.
#' @param nsol Integer. The number of tuning parameters `lambda` at which to
#'   compute Rt.
#' @param delay_distn in the case of a non-gamma delay distribution,
#'   a vector of delay probabilities may be passed here. These will be coerced
#'   to sum to 1, and padded with 0 in the right tail if necessary.
#' @param lambdamin Optional value for the smallest `lambda` to use. This should
#'   be greater than zero.
#' @param lambdamax Optional value for the largest `lambda` to use.
#' @param lambda_min_ratio If neither `lambda` nor `lambdamin` is specified, the
#'   program will generate a lambdamin by lambdamax * lambda_min_ratio.
#'   A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `lambdamin = lambda_min_ratio * lambdamax`.
#'   A very small value will lead to the solution `ct = y`.
#'   This argument has no effect if there is a user-defined `lambda` sequence.
#'
#' @return An object with S3 class `smoothed_mult`. Among the list components:
#' * `y` the observations
#' * `x` a vector of positions at which the counts have been observed.
#' * `convolved` the weighted sum of past `y`.
#' * `ct` the estimated smoothed ratio. This is a matrix with
#'     each column corresponding to one value of `lambda`.
#' * `lambda` the values of `lambda` actually used in the algorithm.
#' * `korder` degree of the estimated piecewise polynomial curve.
#' * `dof` degrees of freedom of the estimated trend filtering problem.
#' * `niter` the required number of iterations for each value of `lambda`.
#' * `convergence` if number of iterations for each value of `lambda` is less
#'     than the maximum number of iterations for the estimation algorithm.
#'
#' @export
#'
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_ct(y)
#' plot(out)
#'
#' out0 <- estimate_ct(y, korder = 0L)
#' plot(out0)
estimate_ct <- function(
    y,
    korder = 3L,
    dist_gamma = c(2.5, 2.5),
    x = 1:n,
    lambda = NULL,
    nsol = 100L,
    delay_distn = NULL,
    lambdamin = NULL,
    lambdamax = NULL,
    lambda_min_ratio = 1e-4,
    maxiter = 1e5,
    init = configure_admm()) {

  arg_is_nonneg_int(korder)
  arg_is_pos_int(nsol, maxiter)
  arg_is_scalar(korder, nsol, lambda_min_ratio)
  arg_is_scalar(lambdamin, lambdamax, allow_null = TRUE)
  arg_is_positive(lambdamin, lambdamax, delay_distn, allow_null = TRUE)
  arg_is_positive(lambda_min_ratio, dist_gamma)
  arg_is_length(2, dist_gamma)
  n <- length(y)

  if (korder + 1 >= n)
    cli_abort("`korder + 1` must be less than observed data length.")

  if (!is.null(delay_distn)) delay_distn <- delay_distn / sum(delay_distn)

  # check that x is a sorted, double vector of length n
  xin <- x
  if (inherits(xin, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)
  arg_is_length(n, x)
  if (is.unsorted(x)) {
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
  }


  convolved <- delay_calculator(y, x - min(x) + 1, dist_gamma, delay_distn)

  if (!inherits(init, "admm_configuration")) {
    cli_abort("`init` must be created with `configure_admm()`.")
  }


  # checks on lambda, lambdamin, lambdamax
  if (is.null(lambda)) lambda <- double(nsol) # prep for create_lambda
  if (is.null(lambdamin)) lambdamin <- -1.0
  if (is.null(lambdamax)) lambdamax <- -1.0
  if (length(lambda) == 0) {
    msg <- "If lambda is not specified,"
    if (lambda_min_ratio >= 1)
      cli_abort("{msg} lambda_min_ratio must be in (0,1).")
    if (lambdamin > 0 && lambdamax > 0 && lambdamin >= lambdamax)
      cli_abort("{msg} lambdamin must be < lambdamax.")
    lambda <- double(nsol)
  }
  if (length(lambda) != nsol) nsol <- length(lambda)
  lambda <- sort(lambda, decreasing = TRUE)

  mod <- estim_path(
    y,
    x,
    convolved,
    korder,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    nsol = nsol,
    rho = init$rho,
    maxiter = maxiter,
    tolerance = init$tolerance,
    lambda_min_ratio = lambda_min_ratio,
    verbose = init$verbose
  )

  structure(
    list(
      y = y,
      x = xin,
      convolved = convolved,
      ct = mod$ct,
      lambda = drop(mod$lambda),
      korder = mod$korder,
      dof = drop(mod$dof),
      niter = drop(mod$niter),
      convergence = (mod$niter < maxiter),
      call = match.call()
    ),
    class = "smoothed_mult"
  )
}


#' ADMM algorithm configuration
#'
#' @param rho Double. An ADMM parameter; coefficient of augmented term in the
#' Lagrangian function.
#' @param tolerance Double. Tolerance of ADMM convergence.
#' @param verbose Integer.
#'
#' @return a list of model parameters with class `admm_configuration`
#'
#' @export
configure_admm <- function(
    rho = -1,
    tolerance = 1e-4,
    verbose = 0) {

  arg_is_scalar(rho, tolerance, verbose)
  arg_is_positive(tolerance)
  arg_is_numeric(rho, tolerance)
  arg_is_nonneg_int(verbose)

  structure(
    list(
      rho = rho,
      tolerance = tolerance,
      verbose = verbose
    ),
    class = "admm_configuration"
  )
}
