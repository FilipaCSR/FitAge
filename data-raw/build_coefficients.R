# =============================================================================
# BUILD THE BUNDLED COEFFICIENT TABLE  ->  data/fitage_coefficients.rda
# =============================================================================
# Combines three kinds of source:
#   1. NHANES joint microdata        -> fit_kdm_coefficients()  (grip, whtr)
#   2. Published grouped summary CSVs -> fit_from_grouped()      (Norway, France)
#   3. Direct k/q/s CSVs              -> used as-is              (reaction time)
#
#   Rscript data-raw/build_coefficients.R
#
# The NHANES step needs internet + the `nhanesA` package; everything else is local.
# =============================================================================

for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) source(f)
shared <- c("bm", "sex", "k", "q", "s", "age_min", "age_max", "source", "notes")

# --- 1. NHANES 2011-2014: grip strength + waist-to-height (joint anchor) ------
build_nhanes <- function() {
  library(nhanesA)
  # Per-hand, per-trial grip readings (kg): 3 trials x 2 hands.
  trial_cols <- c("MGXH1T1", "MGXH1T2", "MGXH1T3",
                  "MGXH2T1", "MGXH2T2", "MGXH2T3")
  fetch <- function(demo, body, grip) {
    d <- nhanes(demo)[, c("SEQN", "RIAGENDR", "RIDAGEYR")]
    b <- nhanes(body)[, c("SEQN", "BMXWAIST", "BMXHT")]
    g <- nhanes(grip)[, c("SEQN", trial_cols)]
    Reduce(function(x, y) merge(x, y, by = "SEQN"), list(d, b, g))
  }
  nh <- rbind(fetch("DEMO_G", "BMX_G", "MGX_G"),
              fetch("DEMO_H", "BMX_H", "MGX_H"))
  nh$age  <- nh$RIDAGEYR
  # nhanesA returns RIAGENDR as a labelled factor ("Male"/"Female"); recode to
  # the 1=male / 2=female convention used by every other source.
  nh$sex  <- ifelse(nh$RIAGENDR %in% c("Male", 1), 1L, 2L)
  # Best SINGLE-HAND grip = max single reading across both hands / all trials.
  # (Matches a home dynamometer "best squeeze" and CLSA dominant-hand norms.)
  trials  <- sapply(nh[trial_cols], as.numeric)
  nh$grip <- apply(trials, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  nh$whtr <- nh$BMXWAIST / nh$BMXHT
  nh <- nh[nh$age >= 18 & nh$age <= 80, ]
  fit_kdm_coefficients(nh, c("grip", "whtr"), source = "NHANES 2011-2014")
}

# --- 2. Grouped summary sources (mean/SD per age band) ------------------------
read_grouped <- function(path, source) {
  df <- read.csv(path, comment.char = "#")
  fit_from_grouped(df, source = source)
}
norway <- read_grouped("data-raw/sources/norwegian_aandstad2016.csv",
                       "Aandstad 2016 (Norway), n=726, ages 20-65")
sitrise <- read_grouped("data-raw/sources/sitrise_araujo2020.csv",
                        "Araujo 2020 SRT (provisional, median-based, 46+)")
# Chair rise: published percentile table -> moments -> grouped fit.
chairrise_raw <- read.csv("data-raw/sources/chairrise_clsa_mayhew2023.csv", comment.char = "#")
chairrise_mom <- percentiles_to_moments(
  chairrise_raw,
  probs = c(p5 = .05, p10 = .10, p20 = .20, p25 = .25, p50 = .50,
            p75 = .75, p80 = .80, p90 = .90, p95 = .95)
)
chairrise <- fit_from_grouped(chairrise_mom, source = "Mayhew 2023 CLSA (Supp App 5)")
# CLSA single-leg balance is censored at 60 s (GAMLSS could not be fit) -> NOT used;
# the Norwegian OLSsum balance above is used instead.
# French max-push-ups: alternative protocol, female slope unreliable -> OFF by default.
# france <- read_grouped("data-raw/sources/pushups_french_nassif2012.csv",
#                        "Nassif 2012 (France)")

# --- 3. Direct k/q/s sources --------------------------------------------------
reaction <- read.csv("data-raw/sources/reaction_ukbiobank_derived.csv",
                     comment.char = "#", stringsAsFactors = FALSE)

# --- assemble -----------------------------------------------------------------
# Try NHANES; fall back gracefully if offline so the rest still builds.
nhanes <- tryCatch(build_nhanes(), error = function(e) {
  warning("NHANES step skipped (", conditionMessage(e),
          "). grip & whtr will be MISSING from the table.")
  NULL
})
# Flag provisional markers (approximate coefficients): floor SRT (median-based,
# Brazilian clinical sample) and reaction time (residual SD ~112 ms is a guess).
norway$provisional    <- FALSE
chairrise$provisional <- FALSE
sitrise$provisional   <- TRUE
reaction$provisional  <- TRUE
sharedp <- c(shared, "provisional")

parts <- list(norway[sharedp], sitrise[sharedp], chairrise[sharedp], reaction[sharedp])
if (!is.null(nhanes)) {
  nhanes$provisional <- FALSE
  parts <- c(list(nhanes[sharedp]), parts)
}

fitage_coefficients <- do.call(combine_coefficients, parts)

dir.create("data", showWarnings = FALSE)
save(fitage_coefficients, file = "data/fitage_coefficients.rda", compress = "xz")
message("Wrote data/fitage_coefficients.rda: ", nrow(fitage_coefficients),
        " rows, markers: ", paste(sort(unique(fitage_coefficients$bm)), collapse = ", "))

# REMAINING REFINEMENTS:
#  - chair_rise: values chart-read from Figure 4; swap in exact numbers from
#    Supplementary Appendix 5 when available.
#  - reaction: residual SD is an approximate placeholder (~112 ms).
#  - grip & whtr: run with internet + nhanesA to populate them.
