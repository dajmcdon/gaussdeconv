
#' Interpolate (or extrapolate) ct estimates to intermediate design points
#'
#' @param object A fitted object produced by `estimate_ct()` or
#'   `cv_estimate_ct()`.
#' @param xout a vector of new positions at which Rt should be produced,
#'   but where counts may not have been observed.
#' @param ... additional arguments passed to methods.
#'
#' @return A vector or matrix of interpolated ct estimates.
#' @export
interpolate_ct <- function(object, xout, ...) {
  if (inherits(xout, "Date")) xout <- as.numeric(xout)
  arg_is_numeric(xout)
  UseMethod("interpolate_ct")
}

#' @export
interpolate_ct.default <- function(object, xout, ...) {
  cli_abort(
    "No interpolation methods exist for objects of class {.cls {class(object)}}."
  )
}
