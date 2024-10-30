#' Parse pedigree-codes.
#'
#' @description
#' Read a pedigree codes definition file and return a duplicated data frame,
#' with entries swapped.
#' @import dplyr
#' @import utils
#' @param path   target pedigree code file name
#' @param header whether or not the file has a header
#' @param sep    whether or not the file has a separatory line
#' @param ...    additional arguments for read.table
#' @export
parse_pedigree_codes <- function(path, header = FALSE, sep = "\t", ...) {
  ped_codes <- utils::read.table(
    file      = path,
    header    = header,
    sep       = sep,
    col.names = c("Relationship", "Ind1", "Ind2", "R"),
    ...
  )
  # --- Duplicate rows, by swapping Ind1 Ind2 values
  ped_codes <- rbind(
    ped_codes,
    dplyr::rename(ped_codes, Ind1 = Ind2, Ind2 = Ind1)
  )
}