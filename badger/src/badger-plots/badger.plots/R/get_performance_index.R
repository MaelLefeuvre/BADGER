#' Estimate Classification Accuracy for a set of BADGER results.
#'
#' @description
#' Compute the confusion matrix and classification performance index of a
#' BADGER run for a given method.
#'
#' - The input data frame must at least contain a 'theoric_r' column.
#'   This column contains the theoretical observed relatedness coefficient for
#'   a pair of individuals. e.g. if the tool estimated this pair as Sibling,
#'   the theoric_r must be 0.5
#'
#' - The input data frame must also contain a 'truth_column' (specifiable as
#'   an argument). This column contains the *expected* theoretical relatedness
#'   coefficient for a given pair.
#'
#' @details
#' Classification performance may be evaluated using two alternate methods:
#'
#' - Ordinal classification index (OCI). See (Cardoso & Sousa, 2011).
#'   DOI: https://doi.org/10.1142/S0218001411009093
#'
#' - Uniform Ordinal Classification Index (UOC). See (Silva, et al. 2018).
#'   DOI: https://doi.org/10.1109/IJCNN.2018.8489327
#'
#' @param data_frame a dataframe. Must contain at least a 'theoric_r'
#'        and 'expected_r' column.
#' @param truth_column column name to use as the expected result
#' @param levels expected levels of the matrix.
#' @param ratio rescale values to balance out class values.
#' @param split logical. When TRUE, compute one OCI per pedigree replicate and
#' return the average index + a confidence interval (computed from the sd)
#' @param conf specify the value of the calculated confidence interval
#' @param method specify which performance index to calculate
#' @param gamma specify the value of the gamma parameter. Only useful when
#'        method = "oci"
#' @param beta specify the value of the beta parameter.
#' @return a named list containing both the generated confusion matrix
#' and the ordinal classification index computed from it.
#' - `confusion_matrix: confusion matrix, created from the provided data frame
#' - `value`: classification performance index, (see parameter `method`)
#' - `error`: confidence interval. see `conf` for the value of this interval.
#' - `method`: character string defining the method used to compute `value`
#'
#' @export
get_performance_index <- function(
  data_frame,
  truth_column     = "true_r",
  levels           = c(0.0, 0.125, 0.25, 0.5, 1.0),
  ratio            = FALSE,
  split            = TRUE,
  conf             = 0.95,
  method           = c("oci", "uoc"),
  gamma            = 1.0,
  beta             = ifelse(method == "uoc", NA, 0.75)
) {
  required_columns <- c(truth_column, "expected_r")
  if (!all(required_columns %in% colnames(data_frame))) {
    msg <- paste(
      "data frame is missing any or all required columns: ",
      toString(required_columns)
    )
    stop(msg)
  }

  observed_r <- factor(data_frame$theoric_r,       levels = levels)
  expected_r <- factor(data_frame[[truth_column]], levels = levels)

  cm <- table(expected_r, observed_r)
  if (ratio) {
    cm <- prop.table(cm, margin = 1L) * (sum(cm) / length(levels))
  }

  # ---- If requested, compute one OCI per pedigree replicate, and return the
  # average OCI + a confidence interval (computed from the sd)
  z     <- badger.plots::.critical_z_value(conf)
  n     <- length(unique(data_frame$rep))
  if (split) {
    values <- unlist(
      dplyr::group_split(data_frame, rep) %>% lapply(FUN = function(x) {
        performance <- get_performance_index(
          data_frame   = x,
          truth_column = truth_column,
          levels       = levels,
          ratio        = ratio,
          split        = FALSE,
          conf         = conf,
          method       = method,
          gamma        = gamma,
          beta         = beta
        )
        performance$value
      })
    )
    value   <- mean(values)
    error <- z * stats::sd(values) / sqrt(n)
  } else {
    # If not, directly compute OCI from the merged dataset, and
    # approximate the confidence interval using normal approximation.
    value   <- switch(method[1L],
      oci = cm_to_oci(cm, gamma = gamma, beta = beta),
      uoc = cm_to_uoc(cm, beta = beta)$value
    )
    error <- z * sqrt((value * (1L - value)) / sum(cm))
  }

  list(confusion_matrix = cm, value = value, error = error, method = method[1L])
}

#' Compute the Ordinal Classification Index of a confusion matrix
#'
#' @description
#' Compute classification performance from a confusion matrix using the ordinal
#' classification index metric.
#'
#' @details
#' This function is adapted from pseudo-code found within the original
#' publication describiing the OCI metric. See (Cardoso & Sousa, 2011).
#' DOI: https://doi.org/10.1142/S0218001411009093
#'
#' @keywords internal
#' @param cm a confusion matrix of k * k classes
#' @param gamma gamma parameter for the OCI formula
#' @param beta factor for the beta parameter. Note that beta will equate
#'        to beta = beta / (N*(k-1)^gamma)
#' @return a floating point value, defining the OCI value of the provided
#'         confusion matrix.
cm_to_oci <- function(cm, gamma = 1.0, beta = 0.75) {
  n    <- sum(cm)
  k    <- length(colnames(cm))
  beta <- beta / (n * (k - 1.0)^gamma)

  helper_m2   <- matrix(0L, k, k, dimnames = list(row.names(cm), row.names(cm)))
  gamma_table <- helper_m2

  for (row in 1L:k) {
    for (col in 1L:k) {
      helper_m2[row, col]   <- cm[row, col] * ((abs(row - col))^gamma)
      gamma_table[row, col] <- gamma
    }
  }

  total_dispersion <- sum(helper_m2)^(1.0 / gamma)
  helper_m1        <- cm / (total_dispersion + n)

  error_matrix <- matrix(0.0, k, k,
    dimnames = list(row.names(cm), row.names(cm))
  )
  error_matrix[1L, 1L] <- 1.0 - helper_m1[1L, 1L] + beta * helper_m2[1L, 1L]

  for (row in 2L:k){
    col <- 1L
    error_matrix[row, col] <- error_matrix[row - 1L, col] -
      helper_m1[row, col] +
      beta * helper_m2[row, col]
  }

  for (col in 2L:k) {
    row <- 1L
    error_matrix[row, col] <- error_matrix[row, col - 1L] -
      helper_m1[row, col] +
      beta * helper_m2[row, col]
  }

  for (col in 2L:k) {
    for (row in 2L:k) {
      costup      <- error_matrix[row - 1L, col]
      costleft    <- error_matrix[row, col - 1L]
      lefttopcost <- error_matrix[row - 1L, col - 1L]
      aux         <- min(c(costup, costleft, lefttopcost))

      error_matrix[row, col] <- (
        aux - helper_m1[row, col] + beta * helper_m2[row, col]
      )
    }
  }

  oci <- error_matrix[k, k]
  return(oci)
}

#' Compute the Uniform Ordinal Classification Index of a confusion matrix
#'
#' @description
#' Compute the Uniform Ordinal Classification Index (UOC) from a provided
#' confusion matrix
#'
#' @details
#' This function is adapted from the original publication describing the UOC
#' metric (Silva et al, 2018) and python code provided through (Albuquerque
#' et al. 2021).
#'
#' - (Silva et al. 2018) A Uniform Performance Index for Ordinal Classification
#'   with Imbalanced Classes," 2018 International Joint Conference on Neural
#'   Networks (IJCNN), Rio de Janeiro, Brazil, 2018, pp. 1-8
#'   https://doi.org/10.1109/IJCNN.2018.8489327
#'
#' - (Albuquerque et al. 2021) Ordinal losses for classification of cervical
#'   cancer risk. PeerJ Computer Science 7:e457
#'   https://doi.org/10.7717/peerj-cs.457
#'
#' @keywords internal
#' @param cm a confusion matrix of k * k classes
#' @param gamma gamma parameter for the OCI formula
#' @param beta factor for the beta parameter. Note that beta will equate
#'        to beta = beta / (N*(k-1)^gamma)
#' @return a floating point value in the range (0, 1), defining the OCI value
#'         of the provided confusion matrix.
cm_to_uoc <- function(cm, beta = NA) {
  k           <- NROW(cm) # Number of classes
  gamma       <- 1.0
  per_class_n <- rowSums(cm)

  # ---- Compute AUOC (Area under curve of beta[0, 1)] if beta is set to NA:
  if (is.na(beta[1L])) {
    return(stats::integrate(lower = 0.0, upper = 1.0, f = function(x) {
      cm_to_uoc(cm, beta = x)
    }))
  }

  # ---- compute total dispersion
  helper_m1 <- matrix(0L, nrow = NROW(cm), ncol = NCOL(cm))
  helper_m2 <- matrix(0L, nrow = NROW(cm), ncol = NCOL(cm))

  for (row in seq_len(k)) {
    for (col in seq_len(k)) {
      n                   <- per_class_n[row]
      helper_m2[row, col] <- (cm[row, col] / (n + 0L^n)) *
        (abs(row - col))^gamma
    }
  }

  total_dispersion <- sum(helper_m2)^(1.0 / gamma)
  for (row in seq_len(k)) {
    for (col in seq_len(k)) {
      n                   <- per_class_n[row]
      helper_m1[row, col] <- cm[row, col] / (n + 0L^n)
    }
  }

  # ---- Compute missing classes and adjust k accordingly.
  k_adjust  <- k - sum(rowSums(cm) == 0L)
  helper_m1 <- helper_m1 / (total_dispersion + k_adjust)

  # ---------------------- #
  sapply(beta, FUN = function(x) {
    beta <- x / k_adjust

    # Create error matrix and fill first entry
    error_mat <- matrix(0L, nrow = NROW(cm), ncol = NCOL(cm))
    error_mat[1L, 1L] <- 1.0 - helper_m1[1L, 1L] + beta * helper_m2[1L, 1L]

    # ---- fill first column
    for (r in 2L:k) {
      error_mat[r, 1L] <- (
        error_mat[r - 1L, 1L] - helper_m1[r, 1L] + beta * helper_m2[r, 1L]
      )
    }

    # ---- fill first row
    for (c in 2L:k) {
      error_mat[1L, c] <- (
        error_mat[1L, c - 1L] - helper_m1[1L, c] + beta * helper_m2[1L, c]
      )
    }

    # ---- fill rest of the matrix
    for (c in 2L:k) {
      for (r in 2L:k) {
        costup          <- error_mat[r - 1L, c]
        costleft        <- error_mat[r, c - 1L]
        costlefttop     <- error_mat[r - 1L, c - 1L]
        aux             <- min(costup, costleft, costlefttop)
        error_mat[r, c] <- aux - helper_m1[r, c] + beta * helper_m2[r, c]
      }
    }

    uoc <- error_mat[NROW(cm), NCOL(cm)]
    return(uoc)
  })
}