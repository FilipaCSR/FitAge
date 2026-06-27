# =============================================================================
# CALIBRATION FROM PUBLISHED SUMMARY TABLES
# =============================================================================
# UK Biobank and CLSA individual data are access-controlled; published norms
# come as grouped means/SDs or as percentiles by age band. These helpers turn
# such tables into the same (k, q, s) coefficient schema as fit_kdm_coefficients(),
# so summary-sourced markers slot in next to the NHANES-fitted ones.
# =============================================================================

#' Fit KDM coefficients from grouped means and SDs
#'
#' Regresses age-band means on the band mid-age (weighted by band n) to get the
#' slope `k` and intercept `q`, and pools the within-band SDs into the residual
#' SD `s`. This recovers the individual-level age regression that KDM needs,
#' assuming band means lie on the age line (the usual normative-table case).
#'
#' @param df Data frame with columns `bm`, `sex`, `age_mid`, `mean`, `sd`, `n`.
#' @param source Provenance string.
#' @return A coefficient table (see [validate_coefficients()]).
#' @export
fit_from_grouped <- function(df, source = "published norms") {
  need <- c("bm", "sex", "age_mid", "mean", "sd", "n")
  stopifnot(all(need %in% names(df)))

  rows <- list()
  for (key in unique(paste(df$bm, df$sex, sep = "\r"))) {
    parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
    sub <- df[df$bm == parts[1] & as.character(df$sex) == parts[2], ]
    if (nrow(sub) < 2) {
      warning("Skipping ", parts[1], " sex=", parts[2],
              ": need >= 2 age bands to fit a slope.")
      next
    }
    fit <- stats::lm(mean ~ age_mid, data = sub, weights = sub$n)
    # pooled within-band SD = sqrt( sum(n_i * sd_i^2) / sum(n_i) )
    s_pooled <- sqrt(sum(sub$n * sub$sd^2) / sum(sub$n))
    rows[[length(rows) + 1]] <- data.frame(
      bm      = parts[1],
      sex     = as.numeric(parts[2]),
      k       = unname(stats::coef(fit)["age_mid"]),
      q       = unname(stats::coef(fit)["(Intercept)"]),
      s       = s_pooled,
      age_min = min(sub$age_mid),
      age_max = max(sub$age_mid),
      source  = source,
      notes   = sprintf("%d age bands, total n=%d", nrow(sub), sum(sub$n)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  validate_coefficients(out)
  out
}

#' Convert percentile rows to per-band mean and SD (normal approximation)
#'
#' For tables that report percentiles (e.g. CLSA: P5/P10/P20/P50/P80/P90/P95)
#' rather than mean/SD. Mean is taken as the median (P50); SD is estimated by
#' averaging `(P_hi - P_lo) / (z_hi - z_lo)` over symmetric percentile pairs.
#' Note: poor for heavily censored measures (e.g. time-capped single-leg balance).
#'
#' @param df Data frame with `bm`, `sex`, `age_mid`, `n`, and one column per
#'   percentile named like `p10`, `p50`, `p90`.
#' @param probs Named numeric vector mapping percentile column name to its
#'   probability, e.g. `c(p10 = 0.10, p50 = 0.50, p90 = 0.90)`.
#' @return Data frame with `bm`, `sex`, `age_mid`, `mean`, `sd`, `n` ready for
#'   [fit_from_grouped()].
#' @export
percentiles_to_moments <- function(df, probs) {
  pcols <- names(probs)
  stopifnot(all(c("bm", "sex", "age_mid", "n") %in% names(df)),
            all(pcols %in% names(df)),
            "p50" %in% pcols)
  z <- stats::qnorm(probs)

  est_sd <- function(row) {
    lo <- pcols[probs < 0.5]
    hi <- pcols[probs > 0.5]
    pairs <- min(length(lo), length(hi))
    lo <- rev(lo)[seq_len(pairs)]; hi <- hi[order(probs[hi])][seq_len(pairs)]
    sds <- mapply(function(a, b) {
      (as.numeric(row[[b]]) - as.numeric(row[[a]])) / (z[b] - z[a])
    }, lo, hi)
    mean(sds)
  }

  data.frame(
    bm      = df$bm,
    sex     = df$sex,
    age_mid = df$age_mid,
    mean    = df[["p50"]],
    sd      = apply(df, 1, est_sd),
    n       = df$n,
    stringsAsFactors = FALSE
  )
}
