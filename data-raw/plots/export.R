#!/usr/bin/env Rscript
# =============================================================================
# EXPORT PLOT-READY VALUES FOR EVERY FITAGE MARKER  ->  one long CSV
# =============================================================================
# For each (marker, sex) in the bundled coefficient table, emits:
#   - "band" rows  : the published age-band points (mean, SD, n, 50% endpoints)
#                    where the source provides them (Norway, France*, CLSA chair
#                    rise, Araujo sit-rise). NHANES-fitted markers (grip, whtr)
#                    and the direct-coefficient reaction marker have no
#                    band-level data, so only the line is emitted.
#   - "fit" rows   : the fitted line on a 0.5-yr age grid over [age_min, age_max],
#                    plus the middle-50% ribbon. Ribbon SD is interpolated
#                    between band SDs when available, else constant `s`.
#
#   Rscript data-raw/plots/export.R
#
# Output (single long-format file):
#   data-raw/plots/plot_data_all.csv
#   data-raw/plots/slope_all.csv      (one row per marker x sex, for legends)
#
# *French push-ups (Nassif 2012) is intentionally OFF in the FitAge calibration
# (alternative protocol, female slope unreliable) and is therefore excluded
# from this export. Edit `include_french` below to add it for comparison.
# =============================================================================

for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) source(f)

include_french <- FALSE
Z50 <- stats::qnorm(0.75)   # ~0.6745 -> middle 50% under normal approximation
AGE_STEP <- 0.5

# ---- marker metadata --------------------------------------------------------
meta <- data.frame(
  marker       = c("grip", "whtr", "sit_reach", "balance_oneleg",
                   "pushups_mod", "sit_rise_floor", "chair_rise", "reaction"),
  marker_label = c("Grip strength",
                   "Waist-to-height ratio",
                   "Sit-and-reach",
                   "One-leg balance (OLSsum)",
                   "Modified push-ups",
                   "Sit-rise floor test",
                   "5-rep chair rise",
                   "Reaction time"),
  units        = c("kg", "ratio", "cm", "s", "reps / 40 s",
                   "score (0-10)", "s", "ms"),
  higher_is    = c("younger", "older", "younger", "younger",
                   "younger", "younger", "older", "older"),
  stringsAsFactors = FALSE
)
if (include_french) {
  meta <- rbind(meta, data.frame(
    marker       = "pushups",
    marker_label = "Push-ups to exhaustion (French)",
    units        = "reps",
    higher_is    = "younger",
    stringsAsFactors = FALSE
  ))
}

# ---- pull every grouped/percentile source as a uniform band table -----------
read_csv0 <- function(p) read.csv(p, comment.char = "#")

bands_norway   <- read_csv0("data-raw/sources/norwegian_kjaer2016.csv")
bands_sitrise  <- read_csv0("data-raw/sources/sitrise_araujo2020.csv")

# CLSA chair rise is published as percentiles -> convert to mean/SD.
chairrise_raw  <- read_csv0("data-raw/sources/chairrise_clsa_mayhew2023.csv")
bands_chair    <- percentiles_to_moments(
  chairrise_raw,
  probs = c(p5 = .05, p10 = .10, p20 = .20, p25 = .25, p50 = .50,
            p75 = .75, p80 = .80, p90 = .90, p95 = .95)
)

bands_all <- rbind(bands_norway, bands_sitrise, bands_chair)
if (include_french) {
  bands_all <- rbind(bands_all,
                     read_csv0("data-raw/sources/pushups_french_nassif2012.csv"))
}

# Direct k/q/s source (reaction): emitted line-only later.
reaction_direct <- read_csv0("data-raw/sources/reaction_ukbiobank_derived.csv")

# Bundled coefficient table = ground truth for k/q/s/age_min/age_max/source.
load("data/fitage_coefficients.rda")
coefs <- fitage_coefficients
stopifnot(all(meta$marker %in% coefs$bm))

# ---- helpers ----------------------------------------------------------------
get_bands <- function(marker, sex) {
  sub <- bands_all[bands_all$bm == marker & bands_all$sex == sex, ]
  if (nrow(sub) == 0) return(NULL)
  sub <- sub[order(sub$age_mid), ]
  data.frame(
    age   = sub$age_mid,
    n     = sub$n,
    mean  = sub$mean,
    sd    = sub$sd,
    p25   = sub$mean - Z50 * sub$sd,
    p75   = sub$mean + Z50 * sub$sd
  )
}

interp_sd <- function(bands, ages, fallback_s) {
  if (is.null(bands) || nrow(bands) < 2) return(rep(fallback_s, length(ages)))
  approx(x = bands$age, y = bands$sd, xout = ages, rule = 2)$y
}

# ---- assemble one long table ------------------------------------------------
out <- list()
slope_rows <- list()

for (i in seq_len(nrow(meta))) {
  m  <- meta$marker[i]
  ml <- meta$marker_label[i]
  un <- meta$units[i]
  hi <- meta$higher_is[i]

  for (sx in c(1, 2)) {
    cf <- coefs[coefs$bm == m & coefs$sex == sx, ]
    if (nrow(cf) == 0) next

    bands  <- get_bands(m, sx)
    n_total <- if (!is.null(bands)) sum(bands$n) else NA_integer_
    n_bands <- if (!is.null(bands)) nrow(bands)  else 0L

    slope_rows[[length(slope_rows) + 1]] <- data.frame(
      marker       = m,
      marker_label = ml,
      units        = un,
      higher_is    = hi,
      sex_code     = sx,
      sex          = ifelse(sx == 1, "Men", "Women"),
      k            = cf$k,
      q            = cf$q,
      s_pooled     = cf$s,
      age_min      = cf$age_min,
      age_max      = cf$age_max,
      n_total      = n_total,
      n_bands      = n_bands,
      provisional  = cf$provisional,
      source       = cf$source,
      stringsAsFactors = FALSE
    )

    common <- list(
      marker       = m,
      marker_label = ml,
      units        = un,
      higher_is    = hi,
      sex_code     = sx,
      sex          = ifelse(sx == 1, "Men", "Women"),
      source       = cf$source,
      provisional  = cf$provisional,
      k            = cf$k,
      q            = cf$q,
      s_pooled     = cf$s,
      age_min      = cf$age_min,
      age_max      = cf$age_max,
      n_total      = n_total,
      n_bands      = n_bands
    )

    # --- band rows (points) -------------------------------------------------
    if (!is.null(bands)) {
      out[[length(out) + 1]] <- data.frame(
        common,
        row_type = "band",
        age      = bands$age,
        n_band   = bands$n,
        mean     = bands$mean,
        sd       = bands$sd,
        p25      = bands$p25,
        p75      = bands$p75,
        fit      = NA_real_,
        lo50     = NA_real_,
        hi50     = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    # --- fit rows (dense line + ribbon) -------------------------------------
    ages <- seq(cf$age_min, cf$age_max, by = AGE_STEP)
    fit  <- cf$q + cf$k * ages
    sdv  <- interp_sd(bands, ages, fallback_s = cf$s)
    out[[length(out) + 1]] <- data.frame(
      common,
      row_type = "fit",
      age      = ages,
      n_band   = NA_integer_,
      mean     = NA_real_,
      sd       = NA_real_,
      p25      = NA_real_,
      p75      = NA_real_,
      fit      = fit,
      lo50     = fit - Z50 * sdv,
      hi50     = fit + Z50 * sdv,
      stringsAsFactors = FALSE
    )
  }
}

plot_data <- do.call(rbind, out)
slope_all <- do.call(rbind, slope_rows)
rownames(plot_data) <- NULL
rownames(slope_all) <- NULL

# stable column order
col_order <- c("marker", "marker_label", "units", "higher_is",
               "sex_code", "sex", "source", "provisional",
               "row_type", "age",
               "n_band", "mean", "sd", "p25", "p75",
               "fit", "lo50", "hi50",
               "k", "q", "s_pooled",
               "age_min", "age_max", "n_total", "n_bands")
plot_data <- plot_data[, col_order]

# round for readable CSV (still enough precision to redraw)
num_round <- function(x, d) ifelse(is.na(x), NA_real_, round(x, d))
plot_data$mean <- num_round(plot_data$mean, 3)
plot_data$sd   <- num_round(plot_data$sd,   4)
plot_data$p25  <- num_round(plot_data$p25,  3)
plot_data$p75  <- num_round(plot_data$p75,  3)
plot_data$fit  <- num_round(plot_data$fit,  4)
plot_data$lo50 <- num_round(plot_data$lo50, 4)
plot_data$hi50 <- num_round(plot_data$hi50, 4)
plot_data$k    <- round(plot_data$k, 6)
plot_data$q    <- round(plot_data$q, 4)
plot_data$s_pooled <- round(plot_data$s_pooled, 4)

slope_all$k <- round(slope_all$k, 6)
slope_all$q <- round(slope_all$q, 4)
slope_all$s_pooled <- round(slope_all$s_pooled, 4)

outdir <- "data-raw/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(plot_data, file.path(outdir, "plot_data_all.csv"), row.names = FALSE)
write.csv(slope_all, file.path(outdir, "slope_all.csv"),     row.names = FALSE)

cat("Wrote:\n  ",
    file.path(outdir, "plot_data_all.csv"), "  (", nrow(plot_data), " rows)\n  ",
    file.path(outdir, "slope_all.csv"),     "  (", nrow(slope_all), " rows)\n\n",
    sep = "")

cat("Slopes:\n")
print(slope_all[, c("marker", "sex", "k", "q", "s_pooled",
                    "age_min", "age_max", "n_total", "n_bands", "provisional")],
      row.names = FALSE)
