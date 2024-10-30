#' GRUPS-rs Ms to r-coefficient.
#'
#' Convert a GRUPS Ms value to r-coefficient, based on simulated distributions.
#' This operation requires to fetch for ".sims" files within the grups-output
#' results folder. Thus, quite the computationally intensive operation...
#'
#' Behavior:
#'  1. Fetch the .sims files of each pair of individuals
#'  2. For each pair:
#'       1. aggregate simulated PWD according to relatedness order
#'       2. order these averages in decreasing order
#'       3. pick the relatedness order with the highest avg. This value
#'          is expected to belong to the 'Unrelated' order
#'       4. use this value as a normalization value for the observed pwd
#'          of the given pair, and compute an r coefficient estimate from it.
#'          i.e.: 2 * (1 - (obs_pwd / norm_value))
#'
#' @importFrom magrittr %>%
#' @param obs_pwd a vector of observed avg pwd (e.g. Corr.Avg.PWD column)
#' @param pair_name a vector of pair names in the form "Ind1-Ind2"
#' @param path path leading to the corresponding read results folder/tarball
#' @param archived specify whether the path is an archived tarball or not.
#' @return a vector of raw r-coefficient estimates.
grups_ms_to_r <- function(obs_pwd, pair_name, path, archived = TRUE) {
  unrelated_ms <- lapply(
    X = read_maybe_tarball(
      path     = path,
      archived = archived,
      pattern  = ".sims$",
      multiple = TRUE,
      sep      = "\t",
      header   = FALSE
    ),
    FUN = function(dataframe) {
      tapply(dataframe$V9, dataframe$V2, FUN = mean) %>% max()
    }
  )

  # Order of extracted tar balls may not be that of the original dataframe
  reorder <- sapply(
    X   = pair_name,
    FUN = function(pair) which(grepl(pair, names(unrelated_ms)))
  )

  if (length(reorder) != length(pair_name)) {
    stop(
      "Error: Multiple matching .sims files when converting ",
      "GRUPS-rs' Ms value to r-coefficient"
    )
  }

  2.0 * (1.0 - (obs_pwd / unlist(unrelated_ms[reorder])))
}

#' Open a GRUPS '.results' file, return a dataframe.
#' @param file Either a .results file or a tarball
#' @param archived Speceify whether or not the provided file is a tarball
#' @param pattern PCRE for .results file. ignored if archive is FALSE
#' @export
load_grups_results <- function(
  file, archived = TRUE, pattern = ".result$"
) {
  if (!archived)
    pattern <- NULL

  read_maybe_tarball(file, archived, pattern, sep = "\t", header = TRUE)
}

#' Format GRUPS-rs results
#'
#' Read from a directory or tarball (see 'archived') and attempt to parse a
#' grups .result file Return a formatted dataframe
#'
#' @details.
#' This will output the results of GRUPS-rs in the form of a data frame,
#' along with several extrapolated variables:
#'
#' - `raw_r_coef`: Contains the raw relatedness coefficient, computed from both
#'   the pairwise mismatch rate of a given pair (`Corr.Avg.Pwd`) and the average
#'   pairwise mismatch rate of Unrelated individuals within corresponding
#'   pedigree simulations of GRUPS-rs. See `badger.plots::grups_ms_to_r` for
#'   more details.
#'
#' - `Ind1`: Contains the label of the first individual involved within the
#'   pair. Extracted from the column `Pair_name`
#'
#' - `Ind2`: Contains the label of the second individual involved within the
#'   pair. Extracted from the column `Pair_name`
#'
#' - `theoric_r`: Contains the expected theoretical r-coefficient of the pair,
#'   given the relationship label provided by GRUPS-rs
#'   (i.e.: `Most_Likely_rel` column)
#'
#' - `true_r`: Contains the truth value for the given pair. This is simply given
#'   out by merging the results of GRUPS-rs with the pedigree-codes dataframe,
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
#' @importFrom tidyr separate_wider_regex
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr mutate filter rename select
#' @importFrom rlang .data
#' @param path Either a directory containing GRUPS results, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param file_name Expected name of the GRUPS .result file
#' @param sample_regex PCRE used to extract and separate sample names
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return a formatted dataframe, containing raw GRUPS-rs results, along with
#' several extrapolated columns (see details.)
#' @export
format_grups_results <- function(
  path              = NULL,
  ped_codes         = NULL,
  archived          = TRUE,
  file_name         = "*.result$",
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
  # the list with  both orders to ensure a match is found.
  relevant_relationships <- c(paste0(ped_codes$Ind1, ped_codes$Ind2))

  # Parse GRUPS '.result' and separate sample names into two separate columns
  results <- load_grups_results(path, archived) %>%
    dplyr::mutate(
      raw_r_coef = grups_ms_to_r(
        Corr.Avg.PWD, Pair_name, path = path, archived = archived
      )
    ) %>%
    tidyr::separate_wider_regex( # split Pair_name into two sep column
      col = .data$Pair_name,
      patterns = c(Ind1 = sample_regex, Ind2 = sample_regex)
    ) %>%
    dplyr::mutate( # Remove delimiter from first match.
      Ind1 = sub("-$", "", .data$Ind1)
    ) %>%
    dplyr::mutate( # Remove unwanted strings from sample names
      Ind1 = gsub(filter_from_names, "", .data$Ind1),
      Ind2 = gsub(filter_from_names, "", .data$Ind2)
    )

  # Filter out unwanted relationships
  subset_results <- results %>% dplyr::filter(
    paste0(.data$Ind1, .data$Ind2) %in% relevant_relationships
  ) %>%
    # Get the 'theoretical' r coefficient from relationships
    dplyr::mutate(
      theoric_r = .rlabel_to_rcoef(.data$Most_Likely_rel, tool = "GRUPS")
    ) %>%
    # Get the 'true' expected relatedness by merging the pedigree codes r_coef
    # values with our data frame.
    merge(ped_codes, by = c("Ind1", "Ind2")) %>%
    dplyr::rename(true_r = "R") %>%
    dplyr::select(-.data$Relationship) %>%
    # Get the r coefficient that we expect from GRUPS.
    # i.e. consider anything greater than the fourth degree
    # as unrelated, as GRUPS does not go beyond this relatedness order
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
  validate_dataset_row_count(subset_results, ped_codes, path)

  subset_results
}