#' badger.plots' version
#'
#' @description
#' Outputs the current version of package `badger.plots`
#'
#' @examples
#' print(badger.plots::version())
#'
#' @export
#' @returns a version tag.
version <- function() {
  utils::packageVersion("badger.plots")
}
