#' Confidence margin to critical z value
#'
#' @description
#' Return a z threshold from a provided confidence interval
#'
#' @keywords internal
#' @param confidence confidence interval
#' @return a z-score, calculated on a reduced-centered gaussian
#' distribution, from the provided confidence interval
#' @export
.critical_z_value <- function(confidence) {
  stats::qnorm((1.0 - confidence) / 2.0, lower.tail = FALSE)
}
