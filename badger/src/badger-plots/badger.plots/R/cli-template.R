#' (Internal) Serialize function arguments
#'
#' @description
#' Serialise the default arguments and parameters of a function as a named list
#'
#' @details
#' This is merely used by `badger.plots::create_template_yaml` to facilitate
#' the creation of template parameters yml configuration files.
#'
#' Note that this function is not recursive, and will not evaluate default
#' values that are themselves lazily-evaluated expressions. In such cases,
#' the default value is simply returned as NA.
#'
#' @keywords internal
#' @param fun a function
#' @return a named list of arguments and default values.
serialize_fn_args <- function(fun) {
  if (typeof(fun) == "closure") {
    fun <- formals(fun)
  }
  fun %>%
    as.list() %>%
    Filter(function(x) !is.symbol(x), .) %>%
    sapply(FUN = function(x) {
      tryCatch(
        expr = {
          if (typeof(x) == "symbol"){
            print(x)
            x
          } else if (typeof(x) == "list language") {
            serialize_fn_args(x)
          } else {
            eval(x)
          }
        },
        error = function(x) return(NULL)
      )
    })
}

#' (Internal) Template badger-plots parameter configuration file.
#'
#' @description
#' Outputs a template yaml file specifying the default command line
#' arguments and parameters required to use the `plot` module of badger-plots'
#' command line argument.
#'
#' @examples
#' param_yaml <- badger.plots::create_template_yaml()
#' print(param_yaml)
#'
#' @import yaml
#' @param path (optional) path to an input `.yml` file. (typically generated
#' through the `make-input` module)
#' @param pedigree_codes (optional) path to a user-created `pedigree-codes.txt`
#' file. See the `badger.plots::parse_pedigree_codes` function for a summary of
#' the required format.
#' @param output_dir (optional) Path specifying the requested output directory
#' of `badger-plots plot` module.
#' @return a YAML-formatted character string
#' @export
create_template_yaml <- function(
  path = NULL, pedigree_codes = NULL, output_dir = NULL
) {
  input_path     <- ifelse(is.null(path), "", path)
  pedigree_codes <- ifelse(is.null(pedigree_codes), "", pedigree_codes)
  output_dir     <- ifelse(is.null(output_dir), "", output_dir)
  yaml           <- list(
    "input"          = input_path,
    "pedigree-codes" = pedigree_codes,
    "output-dir"     = output_dir,
    "reorder"        = vector("character"),
    "rename"         = vector("character"),
    "exclude"        = vector("character")
  )

  yaml[["performance-plot"]]    <- c(
    serialize_fn_args(badger.plots::plot_oci_performance),
    list(width = 1920L, height = 1080L)
  )

  yaml[["accuracy-plot"]]    <- c(
    serialize_fn_args(badger.plots::plot_accuracy),
    list(width = 1120L, height = 1190L)
  )

  yaml <- yaml::as.yaml(yaml)

  header <- paste0(
    "# Package generated with badger.plots version ",
    badger.plots::version(), "\n"
  )
  yaml <- paste0(header, yaml)
  yaml
}

#' (Internal) Argument parser of badger.plots' `template` CLI module.
#'
#' @description
#' Internal wrapper function, which builds and parses the command line argument
#' interface of the `template` module of badger-plots' command line interface
#'
#' @keywords internal
#' @import optparse
#' @param args raw output command line arguments. Typically the output of
#'        the `commandArgs` base function
#' @param ... Additional parameters send to `optparse::parse_args2`
#' @return a hierarchical list of `optparse::OptionParser` command line args
#' Simple command line using optparse
.template_optparse <- function(args = commandArgs(trailingOnly = TRUE), ...) {
  parser <- optparse::OptionParser(
    usage = "%prog [options]",
    add_help_option = TRUE,
    prog = "badger-plots template",
    description = paste(
      "Generate a template yaml file containing all required and optional",
      "parameters for the 'plot' module"
    )
  )
  parser <- optparse::add_option(parser, c("-i", "--input"),
    type    = "character",
    metavar = "path",
    help    = paste(
      "(Optional) Path to an input file. Either a .Rdata containing",
      "previously concatenated results, or a yaml file referencing archived",
      "results for parsing"
    )
  )
  parser <- optparse::add_option(parser, c("-P", "--pedigree-codes"),
    type    = "character",
    metavar = "path",
    dest    = "pedigree_codes",
    help    = paste(
      "(Optional) Path to an input pedigree code definition file"
    )
  )
  parser <- optparse::add_option(parser, c("-o", "--output-dir"),
    type    = "character",
    metavar = "path",
    dest    = "output_dir",
    help    = paste(
      "(Optional) Path to a target output directory where plots and files",
      "should be saved."
    )
  )
  return(optparse::parse_args2(parser, args, ...)$options)
}
