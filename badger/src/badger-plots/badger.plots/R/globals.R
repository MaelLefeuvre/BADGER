# Quick'n'dirty global variable definition. Purely for linting purposes, and
# specifying tidyverse style variables

# format_correctkin_results
utils::globalVariables(c("ID1", "ID2", "corr..kin..coeff."))

# format_grups_results
utils::globalVariables(c("Corr.Avg.PWD", "Pair_name"))

# format_kin_results
utils::globalVariables(c(
  "Pair", "Ind1", "Ind2", "k1", "k2", "true_r", "theoric_r", "expected_r",
  "Relatedness"
))

# format_read_results
utils::globalVariables(c(
  "Ind1", "Ind2", "Rel", "Relationship", "Norm2AlleleDiff",
  "Normalized2AlleleDifference", "true_r", "theoric_r", "expected_r",
  "PairIndividuals"
))

# format_tkgwv2_results
utils::globalVariables(c(
  "Ind1", "Ind2", "Relationship", "HRC", "Relationship.y", "true_r",
  "theoric_r", "expected_r"
))

# parse_pedigree_codes
utils::globalVariables(c("Ind1", "Ind2"))

# serialize_fn_args
utils::globalVariables(".")

# cli-plot-results
utils::globalVariables(c("filename", "width", "height"))

# main
utils::globalVariables(c(
  ".make_input_optparse", ".plot_optparse", ".template_optparse",
  ".plot_all_results", ".rlabel_to_rcoef"
))