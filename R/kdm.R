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
  marker_prec <- sum((k / s)^2)      # precision contributed by the markers
  den <- marker_prec

  corrected <- !is.null(chrono_age) && !is.null(s_ba2)
  prior_prec <- if (corrected) 1 / s_ba2 else 0
  if (corrected) {
    num <- num + chrono_age / s_ba2
    den <- den + prior_prec
  }

  ba <- num / den
  contributions <- num_terms / den   # each marker's share of the estimate, in years
  names(contributions) <- bm

  # Confidence band: total precision `den` -> SE = 1/sqrt(den) (years).
  se <- 1 / sqrt(den)
  ci95 <- c(low = unname(ba - 1.96 * se), high = unname(ba + 1.96 * se))

  # How much of the estimate is driven by markers vs the chronological-age prior.
  marker_information <- marker_prec / den

  # Provisional markers (optional `provisional` column in coefs).
  n_provisional <- if ("provisional" %in% names(coefs))
    sum(as.logical(coefs$provisional[idx]), na.rm = TRUE) else 0L

  stability <- .stability_flag(length(bm), n_provisional, marker_information)

  list(
    functional_age = unname(ba),
    functional_age_advance = if (is.null(chrono_age)) NA_real_ else unname(ba - chrono_age),
    se = unname(se),
    ci95 = ci95,
    marker_information = unname(marker_information),
    stability = stability,
    contributions = contributions,
    markers_used = bm,
    markers_dropped = dropped,
    n_provisional = n_provisional,
    corrected = corrected
  )
}

# Coarse stability flag from marker count, provisional usage, and the share of
# the estimate driven by markers (vs the chronological-age prior).
.stability_flag <- function(n_markers, n_provisional, marker_information) {
  if (n_markers < 3 || marker_information < 0.25) return("low")
  if (n_provisional > 0 || marker_information < 0.45) return("moderate")
  "high"
}

#' Leave-one-out marker sensitivity
#'
#' Recomputes functional age with each marker removed in turn, to show how much
#' each marker moves the score. Large `delta` means the estimate leans heavily on
#' that marker - worth scrutinising for the provisional ones (floor SRT, reaction
#' time). Same arguments as [functional_age()].
#'
#' @inheritParams functional_age
#' @return A data frame (sorted by absolute influence) with columns: `marker`,
#'   `provisional`, `fa_without` (functional age with that marker removed),
#'   `delta` (`fa_without` minus the full-model functional age).
#' @export
leave_one_out <- function(values, coefs, chrono_age = NULL, s_ba2 = NULL,
                          enforce_age_band = TRUE) {
  full <- functional_age(values, coefs, chrono_age, s_ba2, enforce_age_band)
  used <- full$markers_used
  prov_col <- "provisional" %in% names(coefs)

  rows <- lapply(used, function(m) {
    red <- functional_age(values[setdiff(names(values), m)], coefs,
                          chrono_age, s_ba2, enforce_age_band)
    data.frame(
      marker      = m,
      provisional = if (prov_col) isTRUE(as.logical(coefs$provisional[match(m, coefs$bm)])) else NA,
      fa_without  = red$functional_age,
      delta       = red$functional_age - full$functional_age,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-abs(out$delta)), ]
  rownames(out) <- NULL
  out
}

#' Prior variance for the chronological-age correction (`s_ba2`)
#'
#' The corrected Klemera-Doubal estimate (`BA_EC`) is a precision-weighted
#' average of the marker-based estimate and a chronological-age prior;
#' `s_ba2` is the variance of that prior. Smaller `s_ba2` pulls the estimate
#' more strongly toward chronological age, so few or weak markers degrade to
#' "≈ your real age" instead of producing extreme (even negative) values.
#'
#' This is a deliberate modelling choice: a true KDM `s_ba2` is estimated from a
#' cohort with *all* markers measured per person, which public data does not
#' provide for these fitness measures. The default prior SD of 10 years means
#' functional age is expected to fall within roughly ±20 years of chronological
#' age. Because these fitness markers are individually weak age predictors
#' (high noise-to-slope), this prior carries meaningful weight; with more or
#' stronger markers it automatically matters less. Tune `prior_sd_years` to taste.
#'
#' @param prior_sd_years Expected SD (years) of functional age around CA.
#' @return The prior variance `s_ba2` (a scalar), i.e. `prior_sd_years^2`.
#' @export
prior_s_ba2 <- function(prior_sd_years = 10) {
  stopifnot(is.numeric(prior_sd_years), length(prior_sd_years) == 1,
            prior_sd_years > 0)
  prior_sd_years^2
}
