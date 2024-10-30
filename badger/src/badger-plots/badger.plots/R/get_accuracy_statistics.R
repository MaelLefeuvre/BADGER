#' Normalized RMSD and MBE of relatedness coefficients.
#'
#' @description
#' Computes a normalized estimate of the Root Mean Square Deviation and Mean
#' Bias Error of the estimation of relatedness coefficients for a given kinship
#' estimation method.
#'
#' Normalization is applied by dividing the raw RMSE of a given relationship by
#' the “territory” of its theoretical distribution. This is done by calculating
#' the range separating the midrange found between the mean of the given
#' relationship and the adjacent ones. e.g.:
#'
#' Normalization values are as follows:
#'
#' | Relationship  | average expected r | normalization value |
#' | ------------- | ------------------ | ------------------- |
#' | Unrelated     | 0.0                | 0.0625              |
#' | Third degree  | 0.125              | 0.125               |
#' | Second degree | 0.25               | 0.1875              |
#' | First degree  | 0.5                | 0.375               |
#' | Self / Twins  | 1.0                | 0.25                |
#'
#' This normalization step allows comparison both across estimators and
#' relationship types. Partly inspired by the works of (Wang, et al. 2017).
#'
#' @details
#' - The provided dataframe must contains a `true_r` and `raw_r_coef` columns.
#'
#' - Confidence intervals are calculated using the method described within
#'   (Nicholls A., 2014) and (Nicholls A., 2016)
#'
#' @seealso
#' - Bowen Wang, Serge Sverdlov, Elizabeth Thompson, Efficient Estimation of
#'   Realized Kinship from Single Nucleotide Polymorphism Genotypes, Genetics,
#'   Volume 205, Issue 3, 1 March 2017, Pages 1063–1078,
#'   https://doi.org/10.1534/genetics.116.197004
#'
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 1: The calculation of confidence intervals.
#'   J Comput Aided Mol Des 28, 887–918 (2014).
#'   https://doi.org/10.1007/s10822-014-9753-z
#'
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 2: comparing methods.
#'   J Comput Aided Mol Des 30, 103–126 (2016).
#'   https://doi.org/10.1007/s10822-016-9904-5
#'
#' @param data_frame a dataframe of parsed kinship results.
#'        Must contain 'true_r' and 'raw_r_coef' columns.
#' @param conf      confidence interval.
#' @param normalize When true RMSD and MBE estimates are normalized by the
#'        midrange between the previous and next relatedness order
#' @param na.rm Logical indicating whether `NA` values should be stripped prior
#'        to computations.
#' @return a dataframe, containing normalized root mean square deviation
#' and mean bias error estimates.
#' *Columns*:
#' - `true_r`        (numeric) : observed relatedness order class
#' - `norm_value`    (numeric) : normalisation value used for rmsd and mbe
#' - `rmsd`          (numeric) : normalized root mean square error
#' - `rmsd.upper.ci` (numeric) : upper tail length for rmsd
#' - `rmsd.lower.ci` (numeric) : lower tail length for rmsd
#' - `mbe`           (numeric) : normalized mean bias error
#' - `mbe.upper.ci`  (numeric) : upper tail length for mbe
#' - `mbe.lower.ci`  (numeric) : lower tail length for mbe
#'
#' @export
get_accuracy_statistics <- function(
  data_frame = NULL,
  conf       = 0.95,
  normalize  = TRUE,
  na.rm      = FALSE # nolint: object_name_linter.
) {
  output <- data.frame(
    true_r        = numeric(), # observed class
    norm_value    = numeric(), # normalisation value used for rmsd and mbe
    rmsd          = numeric(), # normalized root mean square error
    rmsd.upper.ci = numeric(), # upper tail length for rmsd
    rmsd.lower.ci = numeric(), # lower tail length for rmsd
    mbe           = numeric(), # normalized mean bias error
    mbe.upper.ci  = numeric(), # upper tail length for mbe
    mbe.lower.ci  = numeric()  # lower tail length for mbe
  )

  # --- Loop along the true relatedness coefficient orders.
  classes <- levels(as.factor(data_frame$true_r))
  for (class in classes) {

    # ---- Extract indices corresponding to the current relatedness class
    #      from the dataframe.
    class_indices <- which(data_frame$true_r == class)

    # ---- Compute a normalization value if the user requested it.
    norm_value <- if (normalize) {
      compute_normalization_value(class, classes = classes)
    } else {
      1.0
    }

    # ---- compute mean bias error and root mean square deviation.
    mbe  <- compute_mbe(
      true       = data_frame$true_r[class_indices],
      pred       = data_frame$raw_r_coef[class_indices],
      norm_value = norm_value,
      conf       = conf,
      na.rm      = na.rm
    )

    rmsd <- compute_rmsd(
      true       = data_frame$true_r[class_indices],
      pred       = data_frame$raw_r_coef[class_indices],
      norm_value = norm_value,
      conf       = conf,
      na.rm      = na.rm
    )

    # ---- bind to output dataframe. (@TODO: this is highly inefficient)
    output <- rbind(output,
      data.frame(
        true_r        = class,
        norm_value    = norm_value,
        rmsd          = rmsd$rmsd,
        rmsd.upper.ci = rmsd$upper.ci,
        rmsd.lower.ci = rmsd$lower.ci,
        mbe           = mbe$mbe,
        mbe.upper.ci  = mbe$upper.ci,
        mbe.lower.ci  = mbe$lower.ci
      )
    )
  }
  output
}

#' (Normalized) Root Mean Square Deviation
#'
#' @description
#' Calculate the Root Mean Square Error between a vector of true and observed
#' predictions. Raw output values may optionally be normalized by a provided
#' set value.
#'
#' @details
#' - Confidence intervals are calculated using the method described within
#'   (Nicholls A., 2014) and (Nicholls A., 2016)
#'
#' @seealso
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 1: The calculation of confidence intervals.
#'   J Comput Aided Mol Des 28, 887–918 (2014).
#'   https://doi.org/10.1007/s10822-014-9753-z
#'
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 2: comparing methods.
#'   J Comput Aided Mol Des 30, 103–126 (2016).
#'   https://doi.org/10.1007/s10822-016-9904-5
#'
#' @keywords internal
#' @param true vector of expected relatedness orders
#' @param pred vector of observed relatedness orders
#' @param conf confidence interval (in ratio)
#' @param norm_value value used to normalize the output RMSD.
#' @param na.rm logical value indicating whether NA values should be stripped
#'        before the computation proceeds.
#' @return a named list containing:
#' - `rmsd`     (numeric) : normalized root mean square deviation
#' - `upper.ci` (numeric) : normalized upper tail length.
#' - `lower.ci` (numeric) : normalized lower tail length.
compute_rmsd <- function(
  true,
  pred,
  conf       = 0.95,
  norm_value = 1.0,
  na.rm      = FALSE # nolint: object_name_linter.
) {
  if (length(true) != length(pred))
    stop("compute_mse: length of expected and predicted vectors do not match")

  # Don't account for NA when computing n if na.rm = TRUE
  n <- ifelse(na.rm, sum(!is.na(true)), length(true))

  rmsd <- sqrt(sum((true - pred) ^ 2.0, na.rm = na.rm) / n)

  # ---- Optionally normalize by a user-defined value
  rmsd <- rmsd / norm_value

  # ---- Compute upper and lower confidence intervals.
  z        <- badger.plots::.critical_z_value(conf)
  upper_ci <- rmsd * (sqrt(1.0 + ((z * sqrt(2.0)) / (sqrt(n - 1.0)))) - 1.0)
  lower_ci <- rmsd * (1.0 - sqrt(1.0 - ((z * sqrt(2.0)) / (sqrt(n - 1.0)))))

  list(rmsd = rmsd, upper.ci = upper_ci, lower.ci = lower_ci)

}

#' (Normalized) Mean Bias Error
#'
#' @description
#' Calculate the mean bias error between a vector of true and observed
#' predictions. Raw output values may optionally be normalized by a provided
#' set value.
#'
#' @details
#' - Confidence intervals are calculated using the method described within
#'   (Nicholls A., 2014) and (Nicholls A., 2016)
#'
#' @seealso
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 1: The calculation of confidence intervals.
#'   J Comput Aided Mol Des 28, 887–918 (2014).
#'   https://doi.org/10.1007/s10822-014-9753-z
#'
#' - Nicholls, A. Confidence limits, error bars and method comparison in
#'   molecular modeling. Part 2: comparing methods.
#'   J Comput Aided Mol Des 30, 103–126 (2016).
#'   https://doi.org/10.1007/s10822-016-9904-5
#'
#' @keywords internal
#' @param true vector of expected relatedness orders
#' @param pred vector of observed relatedness orders
#' @param conf confidence interval
#' @param norm_value value used to normalize the output RMSD.
#' @param na.rm logical value indicating whether NA values should be stripped
#'        before the computation proceeds.
#' @return a named list containing:
#' - mbe      (numeric) : normalized root mean square deviation
#' - upper.ci (numeric) : normalized upper tail length.
#' - lower.ci (numeric) : normalized lower tail length.
compute_mbe <- function(
  true,
  pred,
  conf       = 0.95,
  norm_value = 1.0,
  na.rm      = FALSE # nolint: object_name_linter.
) {
  if (length(true) != length(pred))
    stop("compute_mbe: length of expected and predicted vectors do not match")

  # ---- Compute mean bias estimation.
  mbe    <- sum(pred - true, na.rm = na.rm) / length(pred)
  # ---- Optionally normalize by a user-defined value
  mbe    <- mbe / norm_value

  # ---- Compute confidence intervals.
  n      <- length(pred)
  z      <- badger.plots::.critical_z_value(conf)
  mbe_ci <- (z * compute_population_sd(pred, na.rm = na.rm) / sqrt(n))

  list(mbe = mbe, upper.ci = mbe_ci, lower.ci = mbe_ci)
}

#' Population standard deviation
#'
#' @description
#' Compute population standard deviation. This simply swaps out the denominator
#' from (n - 1) to n
#' @keywords internal
#' @param x a numerical vector
#' @param na.rm logical value indicating whether NA values should be stripped
#'        before the computation proceeds.
#' @return (numeric): a population standard deviation estimate.
compute_population_sd <- function(
  x,
  na.rm = FALSE # nolint: object_name_linter.
) {
  sqrt(sum((x - mean(x, na.rm = na.rm)) ^ 2.0, na.rm = na.rm) / length(x))
}

#' Relatedness normalization value calculate.
#'
#' @description
#' Compute a normalization value from a vector of relatedness orders.
#' This value is computed as the midrange between the previous and next
#' relatedness order.e.g. for half-siblings:
#' - half-siblings relatedness: r = 0.25
#' - next relatedness: r =0.5
#' - previous relatedness: r = 0.125
#' - Normalization value = 0.5 - 0.125 = 0.375
#'
#' @details
#' This is mainly used by `badger.plots::get_accuracy_statistics` to
#' normalize the root mean square error and mean bias errror estimates
#' of every tested method and degree of relatedness.
#'
#' @keywords internal
#' @param value (factor)      a vector of relatedness orders
#'                            e.g.: 0, 0.125, 0.25, 0.5, 1
#' @param classes (character) an optional list of factor levels for value
#' @return (numeric) a vector of normalization values
compute_normalization_value <- function(value, classes = levels(value)) {
  # ---- Ensure relatedness levels are bounded between 0 and 1.
  classes <- unique(c(0L, classes, 1L))
  value   <- factor(value, levels = classes)

  # ---- get the indices of our index.
  value_index <- as.numeric(value)

  # ---- Get the lower and upper relatedness order,
  # but clamp values that are lower or greater than 0 / 1.
  lower_value <- as.numeric(
    as.character(classes[clamp(value_index - 1.0, lower = 1.0)])
  )
  upper_value <- as.numeric(
    as.character(classes[clamp(value_index + 1.0, upper = length(classes))])
  )

  # ---- Get the midrange value between the previous and next relatedness order.
  midrange <- (upper_value - lower_value) / 2.0
  midrange
}

#' Clamp values within a given interval.
#'
#' @description
#' Take a vector of numeric values and clamp any value greater/lower than the
#' provided upper and lower thresholds, respectively.
#'
#' @keywords internal
#' @param x     (numeric) a vector of numeric values
#' @param lower (numeric) a lower threshold to clamp values within the vector.
#' @param upper (numeric) an upper threshold to clamp calues within the vector.
#' @param na.rm (boolean) indicate whether NA values should be stripped before
#'        the computation proceeds.
#' @return a clamped vector of numeric values.
clamp <- function(
  x,
  lower = -Inf,
  upper = +Inf,
  na.rm = FALSE # nolint: object_name_linter.
) {
  pmin(pmax(x, lower, na.rm = na.rm), upper, na.rm = na.rm)
}
