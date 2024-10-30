#' Macro and weighted average F-beta index.
#'
#' @description
#' Calculate macro and weighted average F-beta performance index scores for a
#' given BADGER run.
#' @param matrix a confusion matrix
#' @param beta recall factor. a beta of 1 gives out the traditional F1-score.
#'        A beta of 2.0 sets the recall as 2-times more important than accuracy
#' @param average get macro, or weighted average. Allowed values are:
#' - `"macro"`: return the macro average.
#' - `"weighted"`: return the weighted per-class average F1-score.
#' @return a numeric f-beta performance index
#' @export
f_beta_score <- function(matrix = NULL, beta = 1.0, average = NULL) {
  # F_b: (1 + beta^2) * (precision * recall) / ((beta^2)*precision + recall))
  precisions <- compute_precision(matrix)
  recalls    <- compute_recall(matrix)
  f1_scores  <- (1.0 + beta^2.0) *
    (precisions * recalls) / ((beta^2.0) * precisions + recalls)

  if (any(is.nan(f1_scores))) {
    warning("NaN values produced. (divide by zero error ?)")
    f1_scores[is.nan(f1_scores)] <- 0.0
  }

  if (is.null(average)) {
    f1_scores
  } else {
    switch(average,
      "macro"    = sum(f1_scores) / length(f1_scores),
      "weighted" = sum(f1_scores * (rowSums(matrix) / sum(matrix)))
    )
  }
}

#' Precision of a confusion matrix
#'
#' @description
#' returns the average precision of a confusion matrix, i.e. TP/(TP+FP)
#' @param matrix a confusion matrix
#' @return a numeric accuracy score
#' @export
compute_precision <- function(matrix) {
  sapply(seq_len(nrow(matrix)), FUN = function(i) {
    matrix[i, i] / sum(matrix[i, ])
  })
}

#' Recall of a confusion matrix
#'
#' @description
#' Returns the average recall of a confusion matrix. i.e.. TP / (TP + FN)
#' @param matrix a confusion matrix
#' @return a numeric recall score
#' @export
compute_recall <- function(matrix) {
  sapply(seq_len(ncol(matrix)), FUN = function(i) {
    matrix[i, i] / sum(matrix[, i])
  })
}
