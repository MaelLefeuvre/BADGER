
#' HRC To r-coefficient
#'
#' @description
#' convert a tkgwv2 phi coefficient to r coefficient
#' @param hrc a tkgwv2 half relatedness coefficient
#' @return an r-coefficient
tkgwv2_hrc_to_r <- function(hrc) {
  2.0 * hrc
}

#' Parse TKGWV2 results.
#'
#' @description
#' Open a TKGWV2_Results.txt file and return a dataframe
#'
#' @param file     Either a TKGWV2_Results file or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param pattern Regex for results file. ignored if archived = FALSE
#' @export
load_tkgwv2_results <- function(
  file, archived = TRUE, pattern = "TKGWV2_results.txt$"
) {
  if (!archived)
    pattern <- NULL

  read_maybe_tarball(file, archived, pattern, sep = "\t", header = TRUE)
}

#' Format TKGWV2 results.
#'
#' @description
#' Read from a directory or tarball (see 'archived') and attempt to parse
#' a TKGWV2_Results.txt file. Return a dataframe.
#'
#' @details.
#' This will output the results of TKGWV2 in the form of a data frame,
#' along with several extrapolated variables:
#'
#' - `raw_r_coef`: Contains the raw relatedness coefficient, extrapolated
#'   from the value of Phi estimated by the method for a given pair (`HRC`).
#'   See `badger.plots::tkgwv2_hrc_to_r` For more details.
#'
#' - `Ind1`: Contains the label of the first individual involved within the
#'   pair. Extracted from the column `Sample1`
#'
#' - `Ind2`: Contains the label of the second individual involved within the
#'   pair. Extracted from the column `Sample1`
#'
#' - `theoric_r`: Contains the expected theoretical r-coefficient of the pair,
#'   given the relationship label provided by TKGWV2 (i.e.: `Relationship`)
#'
#' - `true_r`: Contains the truth value for the given pair. This is simply
#'   given out by merging the results of TKGWV2 with the pedigree-codes data
#'   frame, provided through the `ped_codes` variable.
#'
#' - `expected_r`: Contains the *'expected'* truth value for the given pair.
#'   This variable is merely intended for relationships which are not expected
#'   to be detected by the method. i.e. TKGWV2 is not expected to detect
#'   relatedness past the second degree. Therefore, third-degree relationships
#'   are *expected* to be classified as Unrelated.
#'
#' - `prediction`: prediction outcome, using the value of `true_r`. See
#'   `badger.plots::r_coef_to_prediction` for a detailled summary of each class
#'
#' - `corrected_prediction`: pediction outcome, using the value of `expected_r`
#'   See `badger.plots::r_coef_to_prediction` for a detailled summary of each
#'   class.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter rename select
#' @param path Either a directory of tkgwv2 results, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param file_name Expected name of the tkgwv2 results file
#' @param sample_regex PCRE pattern used to extract and separate sample names
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return a formatted dataframe, containing raw TKGWV2 results, along with
#' several extrapolated columns (see details.)
#' @export
format_tkgwv2_results <- function(
  path              = NULL,
  ped_codes         = NULL,
  archived          = TRUE,
  file_name         = ".*TKGWV2_Results.txt$",
  sample_regex      = "ped[0-9]{1,2}_[-gGbBiIsS0-9]+",
  filter_from_names = "ped[0-9]{1,2}_",
  ...
) {

  if (length(path) > 1L) {
    stop(
      "Found multiple candidate tkgwv2 results file:\n  - ",
      paste(path, collapse = "\n  - ")
    )
  }

  # Extract relevant relationships from pedigree_codes data frame.
  # Note that the order of each individual might be swapped, so we double
  # the list with both orders to ensure a match is found
  relevant_relationships <- c(paste0(ped_codes$Ind1, ped_codes$Ind2))

  # Parse TKGWV2_Results.txt and format sample names columns
  results <- load_tkgwv2_results(path, archived) %>%
    dplyr::rename(Ind1 = "Sample1", Ind2 = "Sample2") %>% # rename columns
    dplyr::mutate( # remove unwanted pattern using filter_from_names
      Ind1 = gsub(filter_from_names, "", Ind1),
      Ind2 = gsub(filter_from_names, "", Ind2)
    )

  # Filter out unwanted relationships
  subset_results <- results %>% dplyr::filter(
    paste0(Ind1, Ind2) %in% relevant_relationships
  ) %>%
    # Get the 'theoretical' r coefficient from relationships
    dplyr::mutate(
      theoric_r = .rlabel_to_rcoef(Relationship, "TKGWV2")
    ) %>%
    # convert HRC to r coefficient
    dplyr::mutate(
      raw_r_coef = tkgwv2_hrc_to_r(HRC)
    ) %>%
    # Get the 'true' expected relatedness by merging the pedigree codes
    # r coefficient values with our df.
    merge(ped_codes, by = c("Ind1", "Ind2")) %>%
    dplyr::rename(true_r = "R") %>%
    dplyr::select(-Relationship.y) %>%
    # Get the r coefficient that we expect to obtain from tkgwv2
    # i.e. consider anything greater than the second degree as Unrelated
    # as tkgwv2 does no go beyond this relatedness order
    dplyr::mutate(
      expected_r = sapply(true_r, function(x) ifelse(x < 0.25, 0.0, x))
    ) %>%
    # Get a prediction using the *true* r coefficieint as a truth value
    dplyr::mutate(
      prediction = mapply(r_coef_to_prediction, theoric_r, true_r)
    ) %>%
    # Get a prediction using the *expected* r coefficient as a truth value
    dplyr::mutate(
      corrected_prediction = mapply(r_coef_to_prediction, theoric_r, expected_r)
    )

  # Emit a warning if the number of rows is below expectations
  validate_dataset_row_count(subset_results, ped_codes, path)
  subset_results
}
