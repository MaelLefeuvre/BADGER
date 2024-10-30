#' P0 to r-coefficient
#'
#' @description
#' Convert a read normalized P0 to a raw observed r coefficient
#' @keywords internal
#' @param normalized_p0 a READ normalized P0 value
read_norm_p0_to_r <- function(normalized_p0) {
  2.0 * (1.0 - normalized_p0)
}

#' Parse Read Results
#'
#' @description
#' Open a READ_results file. return a dataframe
#'
#' @param file     Either a READ_results file, or a tarball (see 'archived')
#' @param archived Specify whether or not the provided file is a tarball
#' @param pattern  Regex for results file. ignored if archived = FALSE
#' @export
load_read_results <- function(
  file, archived = TRUE, pattern = "READ_results$"
) {
  if (!archived)
    pattern <- NULL

  read_maybe_tarball(file, archived, pattern, sep = "\t", header = TRUE)
}

#' Parse meansP0_AncientDNA_normalized.
#'
#' @description
#' Open a meansP0_AncientDNA_normalized file, return a data frame.
#'
#' @param file Either a meansP0_AncientDNA_normalized, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param pattern Regex for means_P0 file. ignored if archived = FALSE
#' @param v2 logical indicating whether or not the file was generated with
#'        READv2
#' @export
parse_read_means_p0 <- function(
  file, archived = TRUE, pattern = "meansP0_AncientDNA_normalized$", v2 = FALSE
) {
  if (!archived)
    pattern <- NULL
  sep <- ifelse(v2, "\t", " ")
  read_maybe_tarball(file, archived, pattern, sep = sep, header = TRUE)
}

#' Format READ-(v1|v2) Results.
#'
#' @description
#' Read from a directory or tarball (see 'archived') and attempt to parse both
#' a `READ_results` and a `meansP0_AncientDNA_normalized` file. Return a merged
#' dataframe.
#'
#' @details.
#' This will output the results of READ in the form of a data frame,
#' along with several extrapolated variables:
#'
#' - `raw_r_coef`: Contains the raw relatedness coefficient, extrapolated
#'   from normalized P0 estimates of the given pair, and the median of all
#'   estimates. See `badger.plots::read_norm_p0_to_r` For more details.
#'
#' - `Ind1`: Contains the label of the first individual involved within the
#'   pair. Extracted from the column `PairIndividuals`
#'
#' - `Ind2`: Contains the label of the second individual involved within the
#'   pair. Extracted from the column `PairIndividuals`
#'
#' - `theoric_r`: Contains the expected theoretical r-coefficient of the pair,
#'   given the relationship label provided by READ (i.e.: `Relationship` or
#'   `Rel` column, depending on the method's version)
#'
#' - `true_r`: Contains the truth value for the given pair. This is simply given
#'   out by merging the results of READ with the pedigree-codes dataframe,
#'   provided through the `ped_codes` variable.
#'
#' - `expected_r`: Contains the *'expected'* truth value for the given pair.
#'   This variable is merely intended for relationships which are not expected
#'   to be detected by the method. i.e. READ-v1 is not expected to detect
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
#' @importFrom tidyr separate_wider_regex
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter rename select
#' @param path       Either a directory containing read results, or a tarball
#' @param archived   Specify whether or not the provided file is a tarball
#' @param v2         logical specifying whether to expect READv2 results.
#' @param file_names Expected names of the READ results and means_p0 files.
#' @param sample_regex PCRE pattern used to extract and separate sample names
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return a formatted dataframe, containing raw READ-(v1|v2) results, along
#' with several extrapolated columns (see details.)
#' @export
format_read_results <- function(
  path       = NULL,
  ped_codes  = NULL,
  archived   = TRUE,
  v2         = FALSE,
  file_names = list(
    results  = ifelse(v2, "Read_Results.tsv", "READ_results"),
    means_p0 = ifelse(v2,
      "meansP0_AncientDNA_normalized_READv2",
      "meansP0_AncientDNA_normalized"
    )
  ),
  sample_regex      = "ped[0-9]{1,2}_[-gGbBiIsS0-9]+",
  filter_from_names = "ped[0-9]{1,2}_",
  ...
) {
  if (archived) {
    paths <- list(results = path, means_p0 = path)
  } else {
    paths <- list(
      results  = file.path(path, file_names$results),
      means_p0 = file.path(path, file_names$means_p0)
    )
  }


  # ---- Extract relevant relationships from pedigree_codes df.
  # Note that the order of each individual might be swapped, so we double
  # the list with both orders to ensure a match is found.
  relevant_relationships <- c(paste0(ped_codes$Ind1, ped_codes$Ind2))

  # ---- Parse READ_results and means_P0_AncientDNA_normalized and merge them
  # into a single dataframe. Then, split sample names into two separate columns
  merged_results <- merge(
    x  = load_read_results(
      file     = paths$results,
      archived = archived,
      pattern  = paste0(file_names$results, "$")
    ),
    y  = parse_read_means_p0(
      file    = paths$means_p0,
      archived = archived,
      pattern = paste0(file_names$means_p0, "$"),
      v2      = v2
    ),
    by = "PairIndividuals"
  ) %>%
    tidyr::separate_wider_regex( # Split PairIndividuals into two columns
      col = PairIndividuals,
      patterns = c(Ind1 = sample_regex, Ind2 = sample_regex)
    ) %>%
    dplyr::mutate( # remove unwanted strings from sample names (e.g. 'ped4_')
      Ind1 = gsub(filter_from_names, "", Ind1),
      Ind2 = gsub(filter_from_names, "", Ind2)
    )


  # filter out unwanted relationships
  subset_results <- merged_results %>% dplyr::filter(
    paste0(Ind1, Ind2) %in% relevant_relationships
  ) %>%
    # get the 'theoric' r coef. from Relationships
    dplyr::mutate(
      theoric_r = .rlabel_to_rcoef(if (v2) Rel else Relationship, "READ")
    ) %>%
    # convert the normalized P0 to an observed rcoefficient
    dplyr::mutate(
      raw_r_coef = read_norm_p0_to_r(
        if (v2) Norm2AlleleDiff else Normalized2AlleleDifference
      )
    ) %>%
    # Get the 'true' expected relatedness by merging the pedigree codes  r_coef
    # values with our df.
    merge(ped_codes, by = c("Ind1", "Ind2")) %>%
    dplyr::rename(true_r = "R") %>%
    dplyr::select(-ifelse(v2, "Relationship", "Relationship.y")) %>%
    # Get the r_coefficient that we expect to obtain from READ
    # i.e. consider anything greater than the second degree as Unrelated,
    # as READ does not go beyond this relatedness order.
    dplyr::mutate(
      expected_r = sapply(true_r, function(x) ifelse(x < 0.25 && !v2, 0.0, x))
    ) %>%
    # Get a prediction using the *true* r_coefficient as a truth value
    dplyr::mutate(
      prediction = mapply(r_coef_to_prediction, theoric_r, true_r)
    ) %>%
    # Get a prediction using the *expected* r_coefficient as a truth value
    dplyr::mutate(
      corrected_prediction = mapply(r_coef_to_prediction, theoric_r, expected_r)
    )

  # Emit a warning if the number of rows is below the expected.
  validate_dataset_row_count(subset_results, ped_codes, path)
  subset_results
}

#' Format READv2 results.
#'
#' @description
#' Read from a directory or tarball (see 'archived') and attempt to parse both
#' a Read_Results.tsv and a meansP0_AncientDNA_normalized_READv2 file.
#' Return a merged dataframe.
#' @importFrom tidyr separate_wider_regex
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter rename select
#' @param path Either a directory containing read results, or a tarball
#' @param archived Specify whether or not the provided file is a tarball
#' @param file_names Expected names of the READ results and means_p0 files.
#' @param sample_regex PCRE pattern used to extract and separate sample names
#' @param filter_from_names PCRE used to filter unwanted str in sample names
#' @param ped_codes A data frame containing pedigree codes.
#' @param ... Further arguments passed to other methods. Does nothing currently
#' @return A dataframe, containing the raw results of READv2, along with several
#' extrapolated columns. See badger.plots::format_read_results for more details
#' @export
format_readv2_results <- function(
  path       = NULL,
  ped_codes  = NULL,
  archived   = TRUE,
  file_names = list(
    results  = "Read_Results.tsv",
    means_p0 = "meansP0_AncientDNA_normalized_READv2"
  ),
  sample_regex      = "ped[0-9]{1,2}_[-gGbBiIsS0-9]+",
  filter_from_names = "ped[0-9]{1,2}_",
  ...
) {
  format_read_results(
    path              = path,
    archived          = archived,
    file_names        = file_names,
    sample_regex      = sample_regex,
    filter_from_names = filter_from_names,
    ped_codes         = ped_codes,
    v2                = TRUE,
    ...
  )
}
