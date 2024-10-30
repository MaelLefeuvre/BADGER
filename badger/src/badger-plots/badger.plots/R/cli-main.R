#' (Internal) Main command line interface usage.
#'
#' @description
#' Quick and dirty command line interface usage message for `badger-plots`' main
#' command line interface entry point.
#'
#' @details
#' This simply lists the available modules of `badger-plots` and stops the
#' execution of the program.
#'
#' @keywords internal
#' @return NULL
.main_usage <- function() {
  msg <- paste0(
    "badger-plots\n",
    "\n",
    "USAGE: badger-plots <MODULE> [options]\n",
    "\n",
    "MODULES:\n",
    "- make-input: Create a .yml specifying inputs for the 'plot' module\n",
    "- template  : Create a .yml specifying params for the 'plot' module\n",
    "- plot      : Generate accuracy and performance comparison plots\n",
    "- help      : Show this help message and exit\n"
  )
  message(msg)
}

#' Main argument parser of badger-plots' command line interface
#'
#' @description
#' This merely parses the first user-provided argument as a module
#' name and delegates further argument parsing to the corresponding
#' module's own command line parser.
#'
#' @examples
#' raw_args <- c(
#'   "make-input",
#'   "--archive-dir", "/some/badger-archive/directory",
#'   "--subdirs", "0.02X,0.04X,0.06X,0.08X,0.1X,0.2X,0.5X,1X",
#'   "--methods", "READ,READv2,TKGWV2"
#' )
#' args <- badger.plots:::main(raw_args)
#'
#' @import optparse
#' @param args raw output command line arguments. Typically the output of
#' the `commandArgs` base function
#' @return Nature of the value is mainly dependant on which module was selected.
#' See the corresponding function for every module:
#' - "help": `badger.plots:::.main_usage`
#' - "make-input": `badger.plots::generate_input_yaml`
#' - "template": `badger.plots::create_template_yaml`
#' - "plot": `badger.plots::plot_all_results`
#' @export
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  switch(args[1L],
    "help"     = .main_usage(),
    "make-input" = {
      args <- .make_input_optparse(args[-1L])
      make_input_yaml(
        path           = args$archive_dir,
        subdirectories = args$subdirs,
        tools          = args$methods
      ) %>% cat()
    },
    "template" = {
      args <- .template_optparse(args[-1L])
      create_template_yaml(
        path = args$input,
        pedigree_codes = args$pedigree_codes,
        output_dir = args$output_dir
      ) %>% cat()
    },
    "plot"     = {
      args           <- .plot_optparse(args = args[-1L])
      threads        <- args$threads
      data           <- args$yaml$input
      pedigree_codes <- args$yaml$`pedigree-codes`
      output_dir     <- args$yaml$`output-dir`
      exclude        <- c("input", "pedigree-codes", "output-dir")
      optargs        <- args$yaml[!names(args$yaml) %in% exclude]
      .plot_all_results(
        data, pedigree_codes, threads, optargs, output_dir
      )
    },
    .main_usage()
  )
}
