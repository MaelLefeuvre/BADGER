#' Parse BADGER input yaml.
#'
#' @description
#' Read a BADGER input yaml file and return a data frame and return an R object
#' @importFrom yaml read_yaml
#' @param file path leading to an input yaml
#' @export
parse_input_yaml <- function(file) {
  yaml::read_yaml(file)
}