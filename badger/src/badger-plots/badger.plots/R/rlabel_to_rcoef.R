
#' Relationship label to r-coefficient mapper.
#'
#' @keywords internal
#' @description
#' Convert a set of relationship labels to their corresponding expected
#' r-coefficient value.
#' @param labels a list of relationship labels
#' @param tool name of the method from which labels were generated.
#' @return a named list of expected r-coefficients
.rlabel_to_rcoef <- function(
  labels,
  tool = c("correctKin", "KIN", "GRUPS", "READ", "TKGWV2")
) {
  tool <- match.arg(tool)
  r <- .labels2rcoef_mapping(tool)[[1L]][labels]
  if (anyNA(r)) {
    stop(
      "Cannot find a theoric r for ", tool, " relatedness value: ",
      toString(unique(labels[r %>% is.na() %>% which()]))
    )
  }
  r
}

#' rlabel_to_rcoef mapper selector
#'
#' @keywords internal
#' @description
#' Select the appropriate rlabel2rcoef mapper, given the specified method.
#' @param tool name of the kinship estimation method
#' @return a list of relationship labels to r-coefficient mappings.
.labels2rcoef_mapping <- function(tool) {
  mappings <- list(
    correctKin = c(
      "NA-unrelated" = 0.0, # Individuals not found in the .rels.tsv file
      "uncertain"    = 0.0,
      "5th"          = 0.0, # 0.03125,
      "4th"          = 0.0, # 0.0625,
      "3rd"          = 0.125,
      "2nd"          = 0.25,
      "1st"          = 0.5,
      "DUP/MZT"      = 1.0,  # nolint: nonportable_path_linter.
      "JOINT DATA"   = 1.0
    ),
    KIN = c(
      "Unrelated"     = 0.0,
      "Third Degree"  = 0.125,
      "Second Degree" = 0.25,
      "Parent-Child"  = 0.5,
      "Siblings"      = 0.5,
      "Identical"     = 1.0
    ),
    GRUPS = c(
      "Unrelated"     = 0.0,
      "Cousins"       = 0.125,
      "Third Degree"  = 0.125,
      "Second Degree" = 0.25,
      "Half-siblings" = 0.25,
      "First Degree"  = 0.5,
      "Siblings"      = 0.5,
      "Twins"         = 1.0,
      "Self"          = 1.0
    ),
    READ = c(
      "Unrelated"                     = 0.0,
      "Third Degree"                  = 0.125,
      "Second Degree"                 = 0.25,
      "First Degree"                  = 0.5,
      "IdenticalTwins/SameIndividual" = 1.0  # nolint: nonportable_path_linter.
    ),
    TKGWV2 = c(
      "Unrelated"             = 0.0,
      "2nd degree"            = 0.25,
      "1st degree"            = 0.5,
      "Same individual/Twins" = 1.0  # nolint: nonportable_path_linter.
    )
  )
  mappings[tool]
}
