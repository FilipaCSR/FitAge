# =============================================================================
# COEFFICIENT TABLE: SCHEMA, VALIDATION, COMBINATION
# =============================================================================
# The coefficient table is the heart of FitAge. One row per (marker x sex),
# regardless of whether the row came from fitting NHANES or from a published
# normative paper. Columns:
#
#   bm        character  marker id, e.g. "grip", "whtr", "pushups"
#   sex       integer    1 = male, 2 = female (NHANES RIAGENDR convention)
#   k         numeric    slope of marker on age (units of marker per year)
#   q         numeric    intercept of marker on age
#   s         numeric    residual SD of the age regression
#   age_min   numeric    lower age bound the calibration is valid for
#   age_max   numeric    upper age bound
#   source    character  provenance, e.g. "NHANES 2011-2014", "Diaz 2014 CLSA"
#   notes     character  free text (protocol, caveats)
# =============================================================================

.coef_required_cols <- c("bm", "sex", "k", "q", "s")

#' Validate a FitAge coefficient table
#'
#' @param coefs Data frame to check.
#' @return Invisibly returns `coefs` if valid; errors otherwise.
#' @export
validate_coefficients <- function(coefs) {
  if (!is.data.frame(coefs)) stop("`coefs` must be a data frame.")
  missing <- setdiff(.coef_required_cols, names(coefs))
  if (length(missing) > 0) {
    stop("`coefs` is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (any(coefs$s <= 0, na.rm = TRUE)) {
    stop("All residual SD values `s` must be positive.")
  }
  if (!all(coefs$sex %in% c(1, 2))) {
    stop("`sex` must be coded 1 (male) or 2 (female); found: ",
         paste(unique(coefs$sex[!coefs$sex %in% c(1, 2)]), collapse = ", "),
         ". (NHANES RIAGENDR returns 'Male'/'Female' - recode to 1/2.)")
  }
  if (anyNA(coefs[.coef_required_cols])) {
    stop("`coefs` has NA in a required column (bm/sex/k/q/s).")
  }
  if (any(duplicated(coefs[c("bm", "sex")]))) {
    stop("Duplicate (bm, sex) rows in `coefs`; each marker/sex must be unique.")
  }
  invisible(coefs)
}

#' Combine coefficient tables from multiple sources
#'
#' Row-binds calibration tables (e.g. NHANES-fitted markers plus literature
#' markers) into one table and validates the result. Later tables override
#' earlier ones for the same (bm, sex).
#'
#' @param ... One or more coefficient data frames.
#' @return A single validated coefficient table.
#' @export
combine_coefficients <- function(...) {
  tables <- list(...)
  common <- Reduce(intersect, lapply(tables, names))
  combined <- do.call(rbind, lapply(tables, function(t) t[common]))
  # last-write-wins on duplicate (bm, sex)
  key <- paste(combined$bm, combined$sex, sep = "\r")
  combined <- combined[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  rownames(combined) <- NULL
  validate_coefficients(combined)
  combined
}
