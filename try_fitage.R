# =============================================================================
# try_fitage.R - quick interactive test of the FitAge functional-age calculator
# =============================================================================
# Run from the package root:
#     Rscript try_fitage.R
# or  source("try_fitage.R")  in an R session, then call fitage(...) yourself.
#
# Edit the example calls at the bottom, or build your own person and run fitage().
# Supply only the markers you have - the rest are skipped automatically.
# =============================================================================

# --- load the package functions + bundled coefficient table ------------------
for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) source(f)
load("data/fitage_coefficients.rda")

# Marker reference (units / protocol the inputs must match):
#   grip           best SINGLE-HAND grip strength, kg
#   whtr           waist circumference / height (same units), ratio
#   pushups_mod    modified (knee) push-ups completed in 40 s, count
#   balance_oneleg one-leg standing score (OLSsum), seconds
#   sit_reach      sit-and-reach distance, cm
#   chair_rise     time for 5 sit-to-stands, seconds
#   sit_rise_floor floor Sitting-Rising Test composite score 0-10 (46+ only)
#   reaction       reaction time, ms

#' Compute and pretty-print functional age for one person.
#'
#' Uses the chronological-age-corrected KDM (BA_EC) by default, which keeps the
#' estimate sensible with few or weak markers. Tune `prior_sd` (years): smaller
#' = stronger pull toward chronological age. Set prior_sd = Inf for the raw,
#' uncorrected marker-only estimate.
#'
#' @param age   chronological age (years)
#' @param sex   "male"/"female" or 1/2
#' @param ...   named marker values, e.g. grip = 30, chair_rise = 12
#' @param prior_sd prior SD (years) of functional age around chronological age
fitage <- function(age, sex, ..., prior_sd = 10) {
  sex_code <- if (is.character(sex)) {
    if (tolower(sex) %in% c("m", "male", "1")) 1 else 2
  } else as.numeric(sex)

  coefs   <- fitage_coefficients[fitage_coefficients$sex == sex_code, ]
  vals    <- unlist(list(...))
  s_ba2   <- if (is.finite(prior_sd)) prior_s_ba2(prior_sd) else NULL

  res <- functional_age(vals, coefs, chrono_age = age, s_ba2 = s_ba2)
  raw <- functional_age(vals, coefs, chrono_age = age)$functional_age  # marker-only

  cat(strrep("-", 60), "\n")
  cat(sprintf("  %s, age %g\n", ifelse(sex_code == 1, "Male", "Female"), age))
  cat(sprintf("  Functional age : %.1f years   (%+.1f vs actual, %s)\n",
              res$functional_age, res$functional_age_advance,
              ifelse(res$functional_age_advance >= 0, "older", "younger")))
  cat(sprintf("  95%% band       : %.1f - %.1f   (SE %.1f yr)\n",
              res$ci95["low"], res$ci95["high"], res$se))
  cat(sprintf("  Stability      : %s   (%d markers, %d provisional, %.0f%% from markers vs age prior)\n",
              toupper(res$stability), length(res$markers_used),
              res$n_provisional, 100 * res$marker_information))
  cat(sprintf("  [raw marker-only estimate, no CA anchor: %.1f]\n", raw))

  cat("\n  Per-marker contribution (years; high = pushes age up):\n")
  contrib <- sort(res$contributions, decreasing = TRUE)
  for (m in names(contrib)) cat(sprintf("    %-15s %+7.2f\n", m, contrib[m]))
  anchor <- res$functional_age - sum(res$contributions)
  cat(sprintf("    %-15s %+7.2f  (regularization toward actual age)\n",
              "age-anchor", anchor))

  cat("\n  Leave-one-out (score if each marker were removed; * = provisional):\n")
  loo <- leave_one_out(vals, coefs, chrono_age = age, s_ba2 = s_ba2)
  for (i in seq_len(nrow(loo)))
    cat(sprintf("    %-15s %+6.2f yr%s\n", loo$marker[i], loo$delta[i],
                ifelse(isTRUE(loo$provisional[i]), "  *", "")))

  if (length(res$markers_dropped))
    cat("\n  Dropped (age outside calibration):",
        paste(res$markers_dropped, collapse = ", "), "\n")
  cat(strrep("-", 60), "\n\n")
  invisible(res)
}

# =============================================================================
# EXAMPLES - edit these or add your own
# =============================================================================

# A 55-year-old woman, full marker set
fitage(age = 55, sex = "female",
       grip = 28, whtr = 0.52, pushups_mod = 8, balance_oneleg = 45,
       sit_reach = 20, chair_rise = 12.5, sit_rise_floor = 7, reaction = 540)

# A 40-year-old man, only a few markers measured
fitage(age = 40, sex = "male",
       grip = 50, whtr = 0.48, pushups_mod = 16)

# Same man but stronger grip + faster chair rise -> should read younger
fitage(age = 40, sex = "male",
       grip = 60, whtr = 0.45, pushups_mod = 22, chair_rise = 8)
