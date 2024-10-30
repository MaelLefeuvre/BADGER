#' Compare an observed and expected r-coefficient.
#' 
#' @description 
#' Compare an observed r_coefficient with an expected r_coefficient and return
#' a prediction. Possibles values are:
#' - `TP`: True Positive
#' - `TN`: True Negative
#' - `FP`: False Positive
#' - `FN``: False Negative
#' - `UD``: Upper Degree
#' - `LD``: Lower Degree
#' @param obs_rel   an observed r coefficient value
#' @param true_rel  an expected r coefficient value
#' @export
r_coef_to_prediction <- function(obs_rel, true_rel) {
  if (anyNA(c(obs_rel, true_rel))) {
    warning("[r_coef_to_prediction]: Found NA during comparisons")
    return(NA)
  }
  if (obs_rel == 0.0 && true_rel == 0.0)
    return("TN")
  else if (obs_rel == true_rel)
    return("TP")
  else if (obs_rel  == 0.0 && true_rel > 0.0)
    return("FN")
  else if (true_rel == 0.0 && obs_rel  > 0.0)
    return("FP")
  else if (obs_rel > true_rel)
    return("UD")
  else if (obs_rel < true_rel)
    return("LD")
}
