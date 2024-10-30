#' Extract and load a set of files from a given tarball.
#'
#' @description
#' Extract and read all files matching a given pattern from a tarball.
#' Return a list of dataframes.
#'
#' @import archive
#' @param tarball name of the tarball in which the file is located
#' @param pattern PCRE file pattern matcher
#' @param perl    use PCRE regex engine.
#' @param named   should the output list should contain the file path as a name
#' @param ...     Additional arguments for read.table
#' @return a list of dataframes, where each dataframe is a file found within
#' the archive and matching the specified pattern.
#' @export
read_tarballed_files <- function(
  tarball = NULL, pattern = NULL, perl = TRUE, named = TRUE, ...
) {

  # Peek into archive and get the full names of files matching the pattern.
  files          <- archive::archive(tarball)$path
  matching_files <- files[which(grepl(pattern, files, perl = perl))]

  # Read each file as a dataframe using 'read_tarballed_file'
  sapply(
    X         = matching_files,
    FUN       = read_tarballed_file,
    simplify  = FALSE,
    USE.NAMES = named,
    tarball   = tarball,
    ...
  )
}

#' Extract and read a single file from a tarball.
#'
#' @description
#' Extract and read a single file from a tarball. Return a dataframe.
#'
#' @seealso
#'  - `badger.plots::read_tarballed_files` to extract a list of files from an
#'     archive.
#' @import archive
#' @param tarball name of the tarball in which the file is located
#' @param file    complete target filepath within the archive.
#' @param ... Additional arguments for read.table
#' @return a dataframe, containing the contents of the matched archive.
#' @export
read_tarballed_file <- function(tarball, file, ...) {
  utils::read.table(
    archive::archive_read(tarball, file = file),
    ...
  )
}

#' Load a file which may or may not be contained within a tar archive.
#'
#' @description
#' Thin utility wrapper for the function `badger.plots::read_tarballed_files`.
#' This function does makes no assumptions as to whether the provided path
#' *actually* points to an archived file or not.
#
#' This also checks whether multiple files were matched within the archive,
#' which is considered as an error in the specific context of extracting BADGER
#' archive files.
#'
#' @param path Path leading to the results. either a file or a tarball.
#' @param archived optionally specify if the path is compressed
#' @param pattern PCRE file pattern matcher
#' @param multiple specify whether or not the function should expect multiple
#'                 candidate files within the tarball.
#' @param ... specify additional arguments for utils::read.table()
#' @return a single dataframe, containing the contents of the matched file.
#' @export
read_maybe_tarball <- function(
  path,
  archived = TRUE,
  pattern  = NULL,
  multiple = FALSE,
  ...
) {
  # Parse dataframe
  if (archived) {
    dfs <- read_tarballed_files(
      tarball = path,
      pattern = pattern,
      ...
    )

    # Ensure only one file was matched...
    if (!multiple) {
      if (length(dfs) > 1L) {
        msg <- paste(
          "Archived tarball contains multiple files matching the pattern:",
          pattern
        )
        stop(msg)
      }
      dfs <- dfs[[1L]]
    }
  } else {
    if (!is.null(pattern))
      warning("provided file pattern will be ignored, as archived == FALSE")

    dfs <- utils::read.table(path, ...)
  }
  dfs
}
