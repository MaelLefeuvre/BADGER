#' Search archived BADGER results within a directory.
#'
#' @description
#' This function acts as a wrapper to the list.files() function, allowing the
#' user to recursively search for archived BADGER results. These results may
#' either be archived as `.tar.xz` files, or simply kept as simple `.txt`
#' files.
#' @keywords internal
#' @param path root directory where results archives are located
#' @param tool name of the targeted method. Must match the archive's directory
#' @param subdir search within a specific subdirectory of path. the specified
#' subdirectory must be provided relative to the provided `path` argument.
#' @param pattern pattern of the target archive's filename.
#' @return a list of archived files
search_archive_files <- function(
  path,
  tool,
  subdir = "",
  pattern = "*ped[0-9]{1,2}.*(\\.tar\\.xz|\\.txt)$"
) {
  list.files(path, pattern, recursive = TRUE, full.names = TRUE) %>%
    Filter(f = function(x) grepl(paste0(subdir, "/"), x, perl = TRUE)) %>%
    Filter(f = function(x) grepl(paste0(tool, "/"), x, perl = TRUE))
}


#' Create a badger.plots input yaml file.
#'
#' @description
#' Generates a yaml file specifying the inputs required for the badger.plots
#' command line interface.
#'
#' @details
#' The generated yaml file is generally expected to be of depth 2, most likely
#' with the hierarchy `biological condition` > `method` > `[result files]`
#'
#' e.g.:
#' ```
#' 0.02X:
#'   READ:
#'    - /work/badger/results/run-000/results/04-kinship/READ/ped1/ped1.tar.xz
#'    - /work/badger/results/run-000/results/04-kinship/READ/ped2/ped2.tar.xz
#'    - /work/badger/results/run-000/results/04-kinship/READ/ped3/ped3.tar.xz
#'      ...
#'   READv2
#'    - /work/badger/results/run-000/results/04-kinship/READ/ped1/ped1.tar.xz
#'    - /work/badger/results/run-000/results/04-kinship/READ/ped2/ped2.tar.xz
#'      ...
#' 0.04X:
#'   ...
#' ```
#' @examples
#' tools          <- c("READ", "READv2")
#' subdirectories <- c("0.02X", "0.04X")
#' path           <- file.path(Sys.getenv("HOME"), "badger-results")
#'
#' input_yaml <- badger.plots::make_input_yaml(path, subdirectories, tools)
#'
#' @keywords internal
#' @import yaml
#' @param path root directory where all results archives are located
#' @param subdirectories a vector of subdirectories. Each subdirectory must
#'        contain the archived results of a given condition
#'        e.g.: 'CEU - 0.08X - 3% contamination'
#' @param tools a vector of targeted methods. Must match at least one directory
#'              within the archive's tree.
#' @return a YAML-formatted string of characters.
#' @export
make_input_yaml <- function(
  path,
  subdirectories = vector("character"),
  tools          = vector("character")
) {

  sapply(subdirectories, simplify = FALSE, FUN = function(subdir) {
    sapply(X = tools, simplify = FALSE, FUN = search_archive_files,
      path = path, subdir = subdir
    )
  }) %>%
    yaml::as.yaml()
}

#' (Internal) Split a CLI string by a given separator.
#'
#' @description
#' Command line interface utility parser. Splits a user-provided list of
#' arguments separated by character provided by `split`
#'
#' @details
#' This is simply used as a callback for `optparse::add_option`
#'
#' @keywords internal
#' @param short short command line flag
#' @param long long command line flag
#' @param value user-provided value
#' @param object an optparse::OptionParserOption object
#' @param split field-separator character.
#' @return a vector of splitted strings
.parse_charlist <- function(short, long, value, object, split = ",") {
  strsplit(value, split = split)[[1L]]
}

#' (Internal) Argument parser of badger.plots' `make-input` CLI module.
#'
#' @description
#' Internal wrapper function, which builds and parses the command line argument
#' interface of the `make-input` module of badger-plots' command line interface
#'
#' @examples
#' raw_args <- c(
#'   "--archive-dir", "/some/badger-archive/directory",
#'   "--subdirs", "0.02X,0.04X,0.06X,0.08X,0.1X,0.2X,0.5X,1X",
#'   "--methods", "READ,READv2,TKGWV2"
#' )
#' args <- badger.plots:::.make_input_optparse(raw_args)
#'
#' @keywords internal
#' @import optparse
#' @param args raw output command line arguments. Typically the output of
#'        the `commandArgs` base function
#' @param ... Further arguments passed to `optparse::parse_args2`
#' @return a hierarchical list of `optparse::OptionParser` command line args.
.make_input_optparse <- function(
  args = commandArgs(trailingOnly = TRUE), ...
) {

  parser <- optparse::OptionParser(
    usage           = "%prog [options]",
    add_help_option = TRUE,
    prog            = "badger-plots make-input",
    description     = paste(
      "Generate a yaml file targeting target archived results files. The yaml",
      "file is hierarchically ordered according to biological conditions and",
      "targeted methods"
    )
  )
  parser <- optparse::add_option(parser,
    opt_str = c("-d", "--archive-dir"),
    type    = "character",
    metavar = "PATH",
    dest    = "archive_dir",
    help    = paste(
      "(Required) Specify the root directory where all of the requested",
      "archived results subdirectories are located"
    )
  )
  parser <- optparse::add_option(parser,
    opt_str  = c("-s", "--subdirs"),
    type     = "list",
    action   = "callback",
    callback = .parse_charlist,
    metavar  = "SUBDIR,SUBDIR,SUBDIR,...",
    help     = paste(
      "(Required) Specify target archive subdirectories. Each subdirectory",
      "should correspond to a given biological condition (eg. 'CEU - 0.04X')",
      "Note that this argument may accept paths that are either absolute",
      "or relative to the main directory specified through '--archive-dir'.",
      "Regular expressions are also handled."
    )
  )
  parser <- optparse::add_option(parser,
    opt_str  = c("-m", "--methods"),
    type     = "list",
    action   = "callback",
    callback = .parse_charlist,
    metavar  = "TOOLNAME,TOOLNAME,TOOLNAME",
    default  = c("correctKin", "GRUPS", "TKGWV2", "KIN", "READv2", "READ"),
    help     = paste(
      "(Optional) Specify which methods' results should be targeted within",
      "the directory tree"
    )
  )

  args <- optparse::parse_args2(parser, args, ...)$options
  if (is.null(args$archive_dir) || is.null(args$subdirs)) {
    optparse::print_help(parser)
    stop("Missing some or all required arguments")
  }
  args
}
