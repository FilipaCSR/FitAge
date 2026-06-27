# =============================================================================
# KLEMERA-DOUBAL FUNCTIONAL AGE
# =============================================================================
# Combines per-marker age calibrations (k, q, s) into a single functional-age
# estimate. This is the same math TrueAge uses for blood biomarkers, but here
# the coefficients can come from different source datasets per marker.
# =============================================================================

#' Estimate functional age from physical performance measures
#'
#' Implements the Klemera-Doubal Method (KDM). Each marker is assumed to follow
#' a linear age relationship `x = q + k * age` with residual standard deviation
#' `s`, estimated independently (see [fit_kdm_coefficients()] for NHANES-derived
#' markers, or supplied from published norms for the rest).
#'
#' The uncorrected estimate (`BA_E`) is:
#' \deqn{BA_E = \frac{\sum_j (x_j - q_j) k_j / s_j^2}{\sum_j (k_j / s_j)^2}}
#'
#' When both `chrono_age` and `s_ba2` are supplied, the chronological-age
#' correction term is added to numerator and denominator (the corrected KDM,
#' `BA_EC`). `s_ba2` is a property of a *joint* calibration cohort, so it is
#' optional: without it the function returns the unbiased `BA_E`.
#'
#' The sign of `k` encodes direction automatically - grip strength and
#' flexibility decline with age (`k < 0`); waist-to-height ratio and reaction
#' time rise (`k > 0`). No per-marker direction flag is needed.
#'
#' @param values Named numeric vector of measured values for one person. Names
#'   must match `bm` entries in `coefs`. Markers missing from `values` are
#'   skipped (KDM degrades gracefully with fewer markers).
#' @param coefs Data frame of coefficients already filtered to this person's
#'   sex, with columns `bm`, `k`, `q`, `s`. See [validate_coefficients()].
#' @param chrono_age Optional chronological age, used only for the corrected
#'   estimate.
#' @param s_ba2 Optional scalar variance of BA around CA in a calibration
#'   cohort. Required together with `chrono_age` for the corrected estimate.
#' @param enforce_age_band If `TRUE` (default) and `coefs` carries `age_min` /
#'   `age_max` columns and `chrono_age` is given, markers whose calibration does
#'   not cover `chrono_age` are dropped (with a warning) rather than
#'   extrapolated. Set `FALSE` to use every marker regardless of age.
#' @return A list with `functional_age`, `functional_age_advance`
#'   (functional minus chronological, `NA` if `chrono_age` absent),
#'   per-marker `contributions` (years each marker pushes the estimate),
#'   `markers_used`, and `markers_dropped` (out-of-band markers, if any).
#' @export
functional_age <- function(values, coefs, chrono_age = NULL, s_ba2 = NULL,
                           enforce_age_band = TRUE) {
  stopifnot(is.numeric(values), !is.null(names(values)))
  validate_coefficients(coefs)

  bm <- intersect(names(values), coefs$bm)
  if (length(bm) == 0) {
    stop("None of the supplied values match a biomarker in `coefs`.")
  }

  # Drop markers whose calibration age range does not cover the person's age.
  dropped <- character(0)
  if (enforce_age_band && !is.null(chrono_age) &&
      all(c("age_min", "age_max") %in% names(coefs))) {
    bidx <- match(bm, coefs$bm)
    out <- chrono_age < coefs$age_min[bidx] | chrono_age > coefs$age_max[bidx]
    out[is.na(out)] <- FALSE
    if (any(out)) {
      dropped <- bm[out]
      warning("Age ", chrono_age, " is outside the calibration range for: ",
              paste(dropped, collapse = ", "), " - dropped. ",
              "Pass enforce_age_band = FALSE to keep them.")
      bm <- bm[!out]
    }
    if (length(bm) == 0) {
      stop("All matched markers are outside their calibration age range for age ",
           chrono_age, ".")
    }
  }

  idx <- match(bm, coefs$bm)
  k <- coefs$k[idx]
  q <- coefs$q[idx]
  s <- coefs$s[idx]
  x <- values[bm]

  num_terms <- (x - q) * (k / s^2)   # per-marker numerator contribution
  num <- sum(num_terms)
  den <- sum((k / s)^2)

  corrected <- !is.null(chrono_age) && !is.null(s_ba2)
  if (corrected) {
    num <- num + chrono_age / s_ba2
    den <- den + 1 / s_ba2
  }

  ba <- num / den
  contributions <- num_terms / den   # each marker's share of the estimate, in years
  names(contributions) <- bm

  list(
    functional_age = unname(ba),
    functional_age_advance = if (is.null(chrono_age)) NA_real_ else unname(ba - chrono_age),
    contributions = contributions,
    markers_used = bm,
    markers_dropped = dropped,
    corrected = corrected
  )
}
