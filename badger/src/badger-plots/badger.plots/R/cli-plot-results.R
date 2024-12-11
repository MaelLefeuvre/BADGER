#' Deserialize archived BADGER results
#'
#' @description
#' Deserialize a set of archived results of a BADGER run into a set of
#' self contained data frames.
#'
#' @import future
#' @import future.apply
#' @import progressr
#' @importFrom plyr ldply
#'
#' @param yaml a named list containing a hierarchically structured mapping of
#'        BADGER archives. This list is typically generated through the
#'        `make-input` module of the badger-plots command line interface, or
#'        with the `badger.plots:::.make_input_optparse` internal function
#' @param pedigree_codes a dataframe specifying target pedigree relationship
#'        codes. See `badger.plots::parse_pedigree_codes` for more information
#'        regarding the format of this data frame.
#' @param threads specify the number of parallel worker threads.
#' @param method specify which performance metric should be used to estimate
#'        classification performance.
#' @param conf specify which interval to use when estimating CIs
#' @return A named list of hierarchically structured results. Note that the
#' structure of every item mirrors that of the list provided with the `yaml`
#' argument:
#' - `data`: a hierarchical list of data frames, containing the binded results
#'   of every kinship estimation method.
#' - `performance`: a list of `badger.plots::get_performance_index` results.
#' - `statistics`: a list of `badger.plots::get_accuracy_statistics` results.
#' @export
deserialize_results <- function(
  yaml,
  pedigree_codes,
  threads = 1L,
  method = c("uoc", "oci"),
  conf = 0.95
) {
  results     <- list()
  performance <- list()
  statistics  <- list()

  # ---- Set Thread pool
  message("Setting ", threads, " core(s)")
  future::plan(future::multicore, workers = as.numeric(threads))

  for (coverage in names(yaml)) {
    for (tool in names(yaml[[coverage]])) {
      message("Parsing ", coverage, " ", tool)

      # ---- Parse all files
      message("  - Loading all results file...")
      files <- yaml[[coverage]][[tool]]
      progressr::with_progress({
        p <- progressr::progressor(along = seq_len(length(files)))
        dataframes <- future.apply::future_lapply(files, function(file, ...) {
          p(sprintf("file=%s", file))
          .match_method_parser(tool)(
            path      = file,
            ped_codes = pedigree_codes,
            archived  = is_archived(file)
          )
        })
      })

      # ---- Add replicate id
      message("  - Adding replicate id...")
      sapply(seq_along(dataframes), FUN = function(i) dataframes[[i]]$rep <<- i)

      # ---- Combine as a single dataframe
      message("  - Combining results into a single dataframe...")
      results[[coverage]][[tool]] <- plyr::ldply(dataframes)

      # ---- Compute performance classification metrics (UOC / OCI)
      levels <- sort(unique(results[[coverage]][[tool]]$true_r))
      message(
        "  - Computing performance classification metrics ",
        "(method='", method[1L], "')..."
      )
      performance[[coverage]][[tool]] <- get_performance_index(
        results[[coverage]][[tool]], levels = levels, method = method[1L]
      )

      # ---- Compute summary statistics (nRMSD + nMBE)
      message("  - Computing summary statistics...")
      statistics[[coverage]][[tool]] <- get_accuracy_statistics(
        data_frame = results[[coverage]][[tool]],
        conf       = conf,
        normalize  = TRUE,
        na.rm      = TRUE
      )
    }
  }

  # Return all results
  message("Parsing completed.")
  list(data = results, performance = performance, statistics = statistics)
}


#' Parser function mapper function.
#'
#' @description
#' Find the appropriate parser function for a given kinship estimation method.
#'
#' @keywords internal
#' @param tool (character) name of the tool
#' @return a parsing function
.match_method_parser <- function(tool) {
  # Tool 2 parser dictionary
  parser <- list(
    "READ"       = badger.plots::format_read_results,
    "READv2"     = badger.plots::format_readv2_results,
    "KIN"        = badger.plots::format_kin_results,
    "GRUPS"      = badger.plots::format_grups_results,
    "TKGWV2"     = badger.plots::format_tkgwv2_results,
    "correctKin" = badger.plots::format_correctkin_results
  )

  # Use of regex is a hacky way to enable flexible tool naming in the yaml
  idx <- sapply(
    X   = names(parser),
    FUN = function(x) grepl(tool, pattern = paste0(x, "$"))
  )
  parser[[which(idx)]]
}

#' Check if the provided path is a tar archive, i.e. ends with '.tar.gz' or
#' other such extensions.
#' @keywords internal
#' @param file (character) path to a file
#' @return boolean specifying whether this file is expected to be an archive.
is_archived <- function(file) {
  extensions <- c(
    "tar",
    unlist(lapply(c(".t", ".tar."), paste, c("gz", "xz"), sep = ""))
  )
  any(endsWith(file, extensions))
}

#' Plotly html converter.
#'
#' @description
#' Save and export a plotly figure into an htmlfile
#'
#' @details
#' Merely acts as a thin wrapper to `htmlwidgets::saveWidget`.
#'
#' @import htmlwidgets
#' @param fig a `plotly::plot_ly` object
#' @param filename output filename of the exported html (note the `.html`
#' extension is added automatically)
#' @param output_dir output directory location of the html file.
#' @param timestamp when TRUE, adds a timestamp as a prefix to the filename
#' (format: YYYYMMDD-<filename>.svg)
#' @param selfcontained Whether to save the HTML as a single self-contained
#' file, or with external resources and placed in an adjacent directory.
#' @return NULL
#' @export
save_plotly_html <- function(
  fig, filename, output_dir = ".", timestamp = TRUE, selfcontained = TRUE
) {
  prefix   <- ifelse(timestamp, paste0(format(Sys.Date(), "%Y%m%d"), "-"), "")
  filename <- paste0(prefix, filename, ".html")
  path     <- file.path(output_dir, filename)
  htmlwidgets::saveWidget(fig, file = path, selfcontained = selfcontained)
}

#' plotly object to SVG converter.
#'
#' @description
#' Convert and export a plotly figure into an SVG file, using the reticulate
#' package.
#'
#' @details
#' Since plotly requires the use of the reticulate package, this function makes
#' use of conda environments. This conda environment is expected to be created
#' and configured upon installing the badger-plots command line interface, but
#' may also be specified using the `BADGER_PLOTS_CONDA_ENV` environment
#' variable.
#'
#' Note that the BADGER_PLOTS_CONDA_ENV is specified in `$R_HOME/etc/Renviron`
#' configuration file.
#'
#' @import plotly
#' @import reticulate
#' @param fig a `plotly::plot_ly` object
#' @param filename output filename of the exported svg (note the `.svg`
#' extension is added automatically)
#' @param output_dir output directory location of the svg file.
#' @param width width of the svg (in px.)
#' @param height height of the svg (in px.)
#' @param scale scaling factor for the svg width & height.
#' @param timestamp when TRUE, adds a timestamp as a prefix to the filename
#' (format: YYYYMMDD-<filename>.svg)
#' @param ... additional arguments for `plotly::save_image`
#' @return NULL
#' @export
save_plotly_svg <- function(
  fig,
  filename,
  output_dir,
  width     = 1920L,
  height    = 1080L,
  scale     = NULL,
  timestamp = TRUE,
  ...
) {
  conda_env <- Sys.getenv("BADGER_PLOTS_CONDA_ENV")
  if (conda_env != "") {
    condaenv <- Sys.getenv("BADGER_PLOTS_CONDA_ENV")
    conda    <- Sys.getenv("CONDA_EXE")
    reticulate::use_condaenv(condaenv, conda)
  } else {
    message(
      "badger-plots does not appear to have a dedicated conda environment.",
      "Saving plots as static images may not be supported."
    )
  }
  prefix   <- ifelse(timestamp, paste0(format(Sys.Date(), "%Y%m%d"), "-"), "")
  filename <- paste0(prefix, filename, ".svg")
  path     <- file.path(output_dir, filename)
  plotly::save_image(
    fig, path, ..., width = width, height = height, scale = scale
  )
}


#' (Internal) `--yaml` callback function
#'
#' @description
#' Convert a path, provided through the --yaml command line argument
#' and directly parse the contents into a list.
#'
#' @details
#' This is merely used as a callback for an optparse::OptionParserOption
#' object.
#'
#' @keywords internal
#' @import yaml
#' @param short short argument flag
#' @param long long argument flag
#' @param value a user-provided path to a yaml file.
#' @param object an optparse::OptionParserOption instance.
#' @return a deserialized yaml file in the form named list object.
.parse_cli_yaml <- function(short, long, value, object) {
  if (!file.exists(file = value)) {
    stop("Path specified through the --yaml argument does not exist")
  }
  yaml::yaml.load_file(value)
}




#' (Internal) Argument parser of badger.plots' `plot` CLI module.
#'
#' @description
#' Internal wrapper function, which builds and parses the command line argument
#' interface of the `plot` module of badger-plots' command line interface
#'
#' @keywords internal
#' @import optparse
#' @param args raw output command line arguments. Typically the output of
#'        the `commandArgs` base function
#' @param ... Additional parameters send to `optparse::parse_args2`
#' @return a hierarchical list of `optparse::OptionParser` command line args
#' Simple command line using optparse
.plot_optparse <- function(args = commandArgs(trailingOnly = TRUE), ...) {
  parser <- optparse::OptionParser(
    usage = "%prog [options]",
    add_help_option = TRUE,
    prog = "badger.plots plot",
    description = "Some example description"
  )
  parser <- optparse::add_option(parser, c("-y", "--yaml"),
    help     = "Path to an input yaml",
    type     = "list",
    action   = "callback",
    callback = .parse_cli_yaml
  )

  parser <- optparse::add_option(parser, c("-@", "--threads"),
    help = "Number of processing cores",
    type = "integer",
    action = "store",
    default = 1L
  )

  args <- optparse::parse_args2(parser, args, ...)$options
  if (is.null(args$yaml)) {
    stop("--yaml argument is required")
  }

  return(args)
}

#' (Internal) Main function of the `plot` module of badger-plots' CLI.
#'
#' @description
#' Generates multiple plots classification performance, RMSD and MBE
#' plots from a provided input yaml file, and pedigree-code definition file.
#'
#' @keywords internal
#' @param data path leading to a `.yml` file specifying the input data.
#' This yaml file is typically generated using the `make-input` module
#' of the badger-plots' command line interface
#' @param pedigree_codes path leading to a pedigree_codes file. This file must
#' be unheaded, tab-separated and contain the following four fields:
#' - "Relationship": User-defined label of the relationship
#' - "Ind1" : ped-sim identifier of the first individual
#' - "Ind2" : pedsim identifier of the second individual
#' - "R": expected relatedness coefficient shared between Ind1 and Ind2
#' @param threads specify the number of worker threads when parsing results.
#' @param optargs a character list of command line arguments. This would
#'        typically be the output of `badger.plots::.plot_optparse` the raw
#'        output of `base::commandArgs`, or a manually constructed vector of
#'        command line arguments.
#' @param output_dir a character string, specifying the path to a directory
#'        where output files (i.e., html and svg plots, .Rdata backups, etc.)
#'        should be saved.
#' @return a list containing the following items:
#' - results: a list containing all the results used to generate each plot. For
#'   a detailled description of the elements of this list, see the return value
#'   of `badger.plots::deserialize_results`
#' - plots: a list of all generated plots.
.plot_all_results <- function(
  data, pedigree_codes, threads, optargs, output_dir = "."
) {

  output_dir <- ifelse(output_dir %in% c("", NA), ".", output_dir)

  message("Parsing pedigree codes into dataframe...")
  pedigree_codes <- badger.plots::parse_pedigree_codes(pedigree_codes)

  # ---- Prepare output directory and data backup
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  data_backup_path <- file.path(output_dir, "data-backup.Rdata")

  # ---- Parse all results
  if (endsWith(data, ".Rdata")) {
    message("Loading specified .Rdata: ", data)
    load(data)
  } else if (file.exists(data_backup_path)) {
    message("Found data backup: ", data_backup_path)
    message("Loading such data backup instead of the provided input path...")
    load(data_backup_path)
  } else {
    message("Parsing archived results using the specified input yaml...")
    results <- badger.plots::deserialize_results(
      badger.plots::parse_input_yaml(data),
      pedigree_codes,
      threads
    )
    # ---- Save results
    message("Saving parsed results in ", data_backup_path)
    save(results, file = data_backup_path)
  }

  # ---- Reorder, rename and exclude methods if requested
  if (!is.null(optargs$reorder) && length(optargs$reorder) > 0L) {
    o <- optargs$reorder
    for (cov in names(results$data)) {
      results$data[[cov]]               <- results$data[[cov]][o]
      results$performance[[cov]]        <- results$performance[[cov]][o]
      results$statistics[[cov]]         <- results$statistics[[cov]][o]
    }
  }

  if (!is.null(optargs$rename) && length(optargs$rename) > 0L) {

    for (cov in names(results$data)) {
      n <- length(names(results$data[[cov]]))
      names(results$data[[cov]])[1L:n]        <- optargs$rename[1L:n]
      names(results$performance[[cov]])[1L:n] <- optargs$rename[1L:n]
      names(results$statistics[[cov]])[1L:n]  <- optargs$rename[1L:n]
    }
  }

  allow_ragged <- optargs$`ragged-input`
  if (!is.null(allow_ragged) && allow_ragged) {
    keep <- unique(c(unlist(lapply(results$data, FUN=function(x) names(x)))))
  } else {
    keep <- names(results$data[[1L]])
  }

  # ---- Exclude requested entries
  keep <- keep[which(!(keep %in% optargs$exclude))]

  oci_results <- lapply(results$performance, FUN = function(x) x[keep])
  oci_results <- lapply(oci_results, FUN = function(x) x[!is.na(names(x))])
  # ---- Plot overall performance plot

  plotlist <- list()
  if (!is.null(optargs$`performance-plot`)) {
    message("Plotting overall performance plot...")
    oci_fig <- do.call(
      what = badger.plots::plot_oci_performance,
      args = c(
        list(oci_results = oci_results),
        optargs$`performance-plot`
      )
    )


    plotlist <- c(plotlist, list(oci = oci_fig))

    with(optargs$`performance-plot`, {
      save_plotly_html(oci_fig, filename, output_dir)
      save_plotly_svg(oci_fig, filename, output_dir, width, height)
    })
  }

  # ---- Plot overall accuracy plot
  if (!is.null(optargs$`accuracy-plot`)) {
    message("Plotting overall nRMSD plot...")

    rmsd_dfs <- lapply(results$statistics, FUN = function(x) x[keep])
    #rmsd_dfs <- lapply(rmsd_dfs, FUN = function(x) x[!is.na(names(x))])

    rmsd_fig <- do.call(what = badger.plots::plot_accuracy,
      args = c(
        list(rmsd_dfs = rmsd_dfs),
        optargs$`accuracy-plot`
      )
    )
    plotlist <- c(plotlist, list(rmsd = rmsd_fig))

    with(optargs$`accuracy-plot`, {
      save_plotly_html(rmsd_fig, filename, output_dir)
      save_plotly_svg(rmsd_fig, filename, output_dir, width, height)
    })
  }

  message("Done")
  list(results = results, plots = plotlist)
}