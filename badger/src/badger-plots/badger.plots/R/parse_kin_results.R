#' Cotterman to r-coefficient
#'
#' @description
#' convert a set of Cottermann coefficients to a raw r coefficient
#'
#' @keywords internal
#' @param k1 k1 cotterman coefficient: i.e. probability of being IBS1
#' @param k2 k2 cotterman coefficient: i.e. probability of being IBS2
#' @return an r-coefficient
.cottermann_to_r_coef <- function(k1, k2) {
  (k1 / 2.0) + k2
}

#' Open a KIN_results.csv file, return a dataframe.
#' @param file Either a KIN_results.csv file, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param pattern PCRE for KIN_results file. ignored if archived = FALSE
#' @export
load_kin_results <- function(
  file, archived = TRUE, pattern = "KIN_results.csv$"
) {
  if (!archived)
    pattern <- NULL

  read_maybe_tarball(file, archived, pattern, sep = "\t", header = TRUE)
}
#' Format KIN results
#'
#' @description
#' Read from a directory or tarball (see 'archived') and attempt to parse a
#' KIN_results.csv file. Return a formatted dataframe.
#'
#' @details.
#' This will output the results of KIN in the form of a data frame,
#' along with several extrapolated variables:
#'
#' - `raw_r_coef`: Contains the raw relatedness coefficient, extrapolated
#'   from columns `k1` and `k2` See `badger.plots::.cottermann_to_r_coef` for
#'   more details.
#'
#' - `Ind1`: Contains the label of the first individual involved within the
#'   pair. Extracted from the column `Pair`
#'
#' - `Ind2`: Contains the label of the second individual involved within the
#'   pair. Extracted from the column `Pair`
#'
#' - `theoric_r`: Contains the expected theoretical r-coefficient of the pair,
#'   given the relationship label provided by KIN
#'   (i.e.: `Relatedness` column)
#'
#' - `true_r`: Contains the truth value for the given pair. This is simply given
#'   out by merging the results of KIN with the pedigree-codes dataframe,
#'   provided through the `ped_codes` variable.
#'
#' - `expected_r`: Contains the *'expected'* truth value for the given pair.
#'   This variable is merely intended for relationships which are not expected
#'   to be detected by the method.
#'
#' - `prediction`: prediction outcome, using the value of ` true_r`. See
#'   `badger.plots::r_coef_to_prediction` for a detailled summary of each class
#'
#' - `corrected_prediction`: pediction outcome, using the value of `expected_r`
#'   See `badger.plots::r_coef_to_prediction` for a detailled summary of each
#'   class.
#'
#' @importFrom tidyr separate_wider_delim
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter rename select
#' @importFrom stringr str_extract
#' @param path Either a directory containing KIN results, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param file_name Expected name of the KIN_results.csv file.
#' @param sample_regex PCRE used to extract and separate sample names.
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return a formatted dataframe, containing raw KIN results, along with
#' several extrapolated columns (see details.)
#' @export
format_kin_results <- function(
  path              = NULL,
  ped_codes         = NULL,
  archived          = TRUE,
  file_name         = "KIN_results.csv$",
  sample_regex      = "ped[0-9]{1,2}_[-gGbBiIsS0-9]+",
  filter_from_names = "ped[0-9]{1,2}_",
  ...
) {
  if (archived) {
    path <- path
  } else {
    path <- file.path(path, file_name)
  }

  # Extract relevant relationships from pedigree_codes data frame.
  # Note that the order of each individual might be swapped, so we double
  # the list with both orders to ensure a match is found.
  relevant_relationships <- c(paste0(ped_codes$Ind1, ped_codes$Ind2))

  # Parse KIN_results.csv and separate sample names into two separate columns
  results <- load_kin_results(path, archived) %>%
    tidyr::separate_wider_delim( # separate by delimiter
      Pair, names = c("Ind1", "Ind2"), delim = "_._"
    ) %>%
    dplyr::mutate( # extract sample name using pattern
      Ind1 = stringr::str_extract(Ind1, sample_regex),
      Ind2 = stringr::str_extract(Ind2, sample_regex)
    ) %>%
    dplyr::mutate( # remove unwanted pattern using filter_from_names
      Ind1 = gsub(filter_from_names, "", Ind1),
      Ind2 = gsub(filter_from_names, "", Ind2)
    )

  # Filter out unwanted relationships
  subset_results <- results %>% dplyr::filter(
    paste0(Ind1, Ind2) %in% relevant_relationships
  ) %>%
    # Get the 'theoric' r coefficient from Relationships.
    dplyr::mutate(
      theoric_r = .rlabel_to_rcoef(Relatedness, "KIN")
    ) %>%
    # Convert cottermann coefficients to a raw r coefficient
    dplyr::mutate(
      raw_r_coef = .cottermann_to_r_coef(k1, k2)
    ) %>%
    # Get the 'true' expected relatedness by merging the pedigree cods r_coef
    # values with our df.
    merge(ped_codes, by = c("Ind1", "Ind2")) %>%
    dplyr::rename(true_r = "R") %>%
    dplyr::select(-Relationship) %>%
    # Get the r_coefficient that we expect to obtain from KIN
    # i.e. consider anything greater than the third degree as Unrelated,
    # as KIN does not go beyond this relatedness order
    dplyr::mutate(
      expected_r = sapply(true_r, function(x) ifelse(x < 0.125, 0.0, x))
    ) %>%
    # Get a prediction using the *true* r_coefficient as a truth value
    dplyr::mutate(
      prediction = mapply(r_coef_to_prediction, theoric_r, true_r)
    ) %>%
    # Get a prediction using the *expected* r_coefficient as a truth value
    dplyr::mutate(
      corrected_prediction = mapply(r_coef_to_prediction, theoric_r, expected_r)
    )

  # Emit a warning if the number of rows is below expectations
  validate_dataset_row_count(subset_results, ped_codes, path)
  subset_results
}