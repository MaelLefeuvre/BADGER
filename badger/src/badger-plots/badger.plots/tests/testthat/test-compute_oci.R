cm <- function(...) {
  cardinality <- sqrt(length(list(...)))
  t(matrix(nrow = cardinality, ncol = cardinality, c(...)))
}


# See Table 1 of W. Silva 2018: https://doi.org/10.1109/IJCNN.2018.8489327
.uoc_test_matrices <- list(
  "A" = list(
    "expect" = list("0.25" = 0.0, "0.75" = 0.0, "NA" = 0.0),
    "cm"     = cm(
      4L, 0L, 0L, 0L,
      0L, 5L, 0L, 0L,
      0L, 0L, 6L, 0L,
      0L, 0L, 0L, 1L
    )
  ),
  "B" = list(
    "expect" = list("0.25" = 0.46, "0.75" = 0.67, "NA" = 0.56),
    "cm"     = cm(
      0L, 4L, 0L, 0L,
      0L, 0L, 6L, 0L,
      0L, 0L, 5L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "C" = list(
    "expect" = list("0.25" = 0.62, "0.75" = 0.71, "NA" = 0.65),
    "cm"     = cm(
      0L, 0L, 4L, 0L,
      0L, 0L, 6L, 0L,
      0L, 0L, 5L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "D" = list(
    "expect" = list("0.25" = 0.56, "0.75" = 0.67, "NA" = 0.61),
    "cm"     = cm(
      0L, 4L, 0L, 0L,
      6L, 0L, 0L, 0L,
      0L, 0L, 5L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "E" = list(
    "expect" = list("0.25" = 0.68, "0.75" = 0.80, "NA" = 0.74),
    "cm"     = cm(
      0L, 4L, 0L, 0L,
      6L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "F" = list(
    "expect" = list("0.25" = 0.56, "0.75" = 0.67, "NA" = 0.61),
    "cm"     = cm(
      0L, 40L, 0L, 0L,
      6L,  0L, 0L, 0L,
      0L,  0L, 5L, 0L,
      0L,  0L, 0L, 3L
    )
  )
)

# ---- See Table 1 of Cardoso & Sousa, 2011
# doi: https://doi.org/10.1142/S0218001411009093
.oci_test_matrices <- list(
  "A" = list(
    "expect" = list("0.25" = 0.00, "0.75" = 0.00),
    "cm"     = cm(
      4L, 0L, 0L, 0L,
      0L, 6L, 0L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "B" = list(
    "expect" = list("0.25" = 0.50, "0.75" = 0.63),
    "cm"     = cm(
      0L, 4L, 0L, 0L,
      0L, 0L, 6L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "C" = list(
    "expect" = list("0.25" = 0.61, "0.75" = 0.78),
    "cm"     = cm(
      0L, 0L, 4L, 0L,
      0L, 0L, 6L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 3L
    )
  ),
  "D" = list(
    "expect" = list("0.25" = 0.65, "0.75" = 0.72),
    "cm"     = cm(
      0L, 4L, 0L, 0L,
      6L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 3L
    )
  )
)

# ---- test badger.plots::cm_to_uoc using .uoc_test_matrices
for (name in names(.uoc_test_matrices)) {
  test_name <- paste("get_performance_index::cm_to_uoc::cm", name, sep = "-")
  expect    <- .uoc_test_matrices[[name]]$expect
  cm        <- .uoc_test_matrices[[name]]$cm
  testthat::test_that(test_name, {
    sapply(names(expect), simplify = FALSE, function(beta) {
      want <- expect[[beta]]
      beta <- ifelse(beta == "NA", NA, as.numeric(beta))
      got  <- badger.plots:::cm_to_uoc(cm, beta)[[1L]]
      testthat::expect_equal(got, want, tolerance = 1e-02)
    })
  })
}

# ---- test badger.plots::cm_to_oci using .oci_test_matrices
for (name in names(.oci_test_matrices)) {
  test_name <- paste("get_performance_index::cm_to_oci::cm", name, sep = "-")
  expect    <- .oci_test_matrices[[name]]$expect
  cm        <- as.table(.oci_test_matrices[[name]]$cm)
  testthat::test_that(test_name, {
    sapply(names(expect), simplify = FALSE, function(beta) {
      want <- expect[[beta]]
      beta <- as.numeric(beta)
      got  <- badger.plots:::cm_to_oci(cm, beta = beta)[[1L]]
      testthat::expect_equal(got, want, tolerance = 1e-02)
    })
  })
}
