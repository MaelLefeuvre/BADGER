#' Parse correctKin results
#'
#' @description
#' Open a CorrectKin `.rels.tsv` file, return a dataframe.
#'
#' @param file Either a .results file or a tarball
#' @param archived Speceify whether or not the provided file is a tarball
#' @param pattern PCRE for .results file. ignored if archive is FALSE
#' @export
load_correctkin_results <- function(
  file, archived = TRUE, pattern = ".rels.tsv$"
) {
  if (!archived)
    pattern <- NULL

  read_maybe_tarball(file, archived, pattern, sep = "\t", header = TRUE)
}

#' Format correctKin results
#'
#' @description
#' Read from a directory or tarball (see 'archived'), attempt to parse a
#' correctKin .rels.tsv file, and format these results (see details).
#'
#' @details.
#' This will output the results of correctKin in the form of a data frame,
#' along with several extrapolated variables:
#' - `theoric_r`: Contains the expected theoretical r-coefficient of the pair,
#'   given the relationship label provided by correectKin
#'   (i.e.: `est..relatedness` column)
#' - `true_r`: Contains the truth value for the given pair. This is simply given
#'   out by merging the results of correctKin with the pedigree-codes dataframe
#'   provided through the `ped_codes` variable.
#' - `expected_r`: Contains the *'expected'* truth value for the given pair.
#'   This variable is merely intended for relationships which are not expected
#'   to be detected by the method. i.e. correctKin cannot detect relationship
#'   past the 4th-degree; therefore, 5th-degree relationship are *expected* to
#'   be classified out as *Unrelated*
#' - `prediction`: prediction outcome, using the value of ` true_r`. See
#'   `badger.plots::r_coef_to_prediction` for a detailled summary of each class
#' - `corrected_prediction`: pediction outcome, using the value of `expected_r`
#'   See `badger.plots::r_coef_to_prediction` for a detailled summary of each
#'   class.
#'
#' @importFrom tidyr separate_wider_regex
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter rename select
#' @importFrom rlang .data
#' @param path Either a directory containing GRUPS results, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param file_name Expected name of the GRUPS .result file
#' @param sample_regex PCRE used to extract and separate sample names
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return a formatted dataframe, containing raw correctKin results, along with
#' several extrapolated columns (see details.)
#' @export
format_correctkin_results <- function(
  path              = NULL,
  ped_codes         = NULL,
  archived          = TRUE,
  file_name         = "*.rels.tsv$",
  sample_regex      = "ped[0-9]{1,2}_[-gGbBiIsS0-9]+",
  filter_from_names = "ped[0-9]{1,2}_",
  ...
) {
  if (archived) {
    path <- path
  } else {
    path <- file.path(path, file_name)
  }

  # Parse correctKIN '.result' and rename columns.
  results <- load_correctkin_results(path, archived) %>%
    dplyr::rename(Ind1 = ID1, Ind2 = ID2)             # Rename Ind1 and Ind2

  # filter ped[0-9]+ from names.
  results <- results %>% dplyr::mutate(
    Ind1 = gsub(filter_from_names, "", .data$Ind1),
    Ind2 = gsub(filter_from_names, "", .data$Ind2)
  )

  # Set raw_r_coef as twice the corr.kin.coeff.
  results <- results %>% dplyr::mutate(raw_r_coef = 2.0 * corr..kin..coeff.)

  # Left-join pedigree codes. (this filters out unwanted relationships, while
  # ensuring non related individuals are still kept within the dataframe)
  rel_pairs <- paste(results$Ind1, results$Ind2)
  for (i in 1L:(NROW(ped_codes) / 2L)) {
    ind1 <- ped_codes[i, "Ind1"]
    ind2 <- ped_codes[i, "Ind2"]

    pair_found <- function(a, b) paste(a, b) %in% rel_pairs

    if (!pair_found(ind1, ind2) && !pair_found(ind2, ind1)) {
      results <- dplyr::bind_rows(
        results,
        data.frame(Ind1 = ind1, Ind2 = ind2, est..relatedness = "NA-unrelated",
          stringsAsFactors = FALSE
        )
      )
    }

  }

  # Get the 'theoretical' r coefficient from relationships
  results <- results %>% dplyr::mutate(
    theoric_r = .rlabel_to_rcoef(.data$est..relatedness, tool = "correctKin")
  ) %>%
    # Get the 'true' expected relatedness by merging the pedigree codes r_coef
    # values with our data frame.
    merge(ped_codes[, c("Ind1", "Ind2", "R")], by = c("Ind1", "Ind2")) %>%
    dplyr::rename(true_r = "R") %>%
    # Get the r coefficient that we expect from correctKin
    # i.e. consider anything greater than the 4th degree as Unr., as correctKin
    # does not go beyond this tie.
    dplyr::mutate(
      expected_r = sapply(.data$true_r, function(x) ifelse(x < 0.0625, 0.0, x))
    ) %>%
    # Get a prediction using the *true* r coefficient as a truth value
    dplyr::mutate(
      prediction = mapply(r_coef_to_prediction, .data$theoric_r, .data$true_r)
    ) %>%
    # Get a prediction using the *expected* r coefficient as a truth value.
    dplyr::mutate(
      corrected_prediction = mapply(
        r_coef_to_prediction, .data$theoric_r, .data$expected_r
      )
    )

  # Emit a warning if the number of rows is below expectations
  validate_dataset_row_count(results, ped_codes, path)
  return(results)
}