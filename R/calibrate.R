# =============================================================================
# CALIBRATION FROM A JOINT DATASET (e.g. NHANES)
# =============================================================================
# For markers measured on the same individuals (grip strength, waist-to-height
# ratio in NHANES 2011-2014) we estimate k, q, s by regressing each marker on
# age, separately by sex. The output is a coefficient table in the standard
# FitAge schema so it can be combined with literature-sourced markers.
# =============================================================================

#' Fit KDM coefficients (k, q, s) from a joint dataset
#'
#' Regresses each marker on age within each sex stratum and returns the slope
#' (`k`), intercept (`q`) and residual SD (`s`) - the three numbers KDM needs.
#'
#' @param data Data frame with one row per individual. Must contain `age`, a
#'   sex column, and one column per marker in `biomarkers`.
#' @param biomarkers Character vector of marker column names to calibrate.
#' @param age Name of the age column. Default `"age"`.
#' @param sex Name of the sex column (values 1 = male, 2 = female).
#'   Default `"sex"`.
#' @param source Provenance string written into the `source` column.
#' @return A coefficient table (see [validate_coefficients()]).
#' @export
fit_kdm_coefficients <- function(data, biomarkers, age = "age", sex = "sex",
                                 source = "NHANES 2011-2014") {
  stopifnot(all(c(age, sex, biomarkers) %in% names(data)))

  rows <- list()
  for (s_val in sort(unique(data[[sex]]))) {
    sub <- data[data[[sex]] == s_val, , drop = FALSE]
    for (bm in biomarkers) {
      df <- data.frame(y = sub[[bm]], age = sub[[age]])
      df <- df[stats::complete.cases(df), , drop = FALSE]
      if (nrow(df) < 10) next
      fit <- stats::lm(y ~ age, data = df)
      rows[[length(rows) + 1]] <- data.frame(
        bm      = bm,
        sex     = s_val,
        k       = unname(stats::coef(fit)["age"]),
        q       = unname(stats::coef(fit)["(Intercept)"]),
        s       = stats::sigma(fit),
        age_min = min(df$age),
        age_max = max(df$age),
        source  = source,
        notes   = sprintf("n=%d", nrow(df)),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  validate_coefficients(out)
  out
}
