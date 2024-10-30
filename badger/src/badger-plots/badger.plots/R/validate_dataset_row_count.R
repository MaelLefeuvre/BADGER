#' Compare the number of rows of a results dataset with that of our pedigree
#' codes table.
#' We are doubling the entries within our pedigree code table to account for
#' switched samples, so we expect that the number of rows in our table is half
#' of our pedigree codes data frame
#' @param results_df a kinship results dataframe
#' @param ped_codes Pedigree codes dataframe
#' @param path original filepath. Only used when emitting a warning.
#' @export
validate_dataset_row_count <- function(results_df, ped_codes, path) {
  # Check if we have the correct number of rows.
  expect <- NROW(ped_codes) / 2.0
  got    <- NROW(results_df)
  if (got < expect) {
    msg <- paste(
      "Some entries seem to be missing in:", path,
      "Expected", expect, "entries after formatting. Got", got
    )
    warning(msg)
  }

}