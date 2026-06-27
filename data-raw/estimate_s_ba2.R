# =============================================================================
# Can s_ba2 be ESTIMATED from data (rather than chosen as a prior)?
# =============================================================================
# The canonical Klemera-Doubal estimator is
#     s_BA^2 = mean_i[(BA_E,i - CA_i)^2] - 1 / Sum_j (k_j/s_j)^2
# which needs a cohort with EVERY marker measured per person. FitAge has none.
#
# This script runs the estimator on the only joint subset available - NHANES
# grip + waist-to-height (2 markers per person) - to show it is degenerate:
# the 2-marker estimation noise (~31-33 yr SD) is larger than the observed
# variance of BA_E around CA, so the estimate comes out NEGATIVE (men ~ -96,
# women ~ -141). Conclusion: s_ba2 cannot be estimated here and is set as a
# deliberate regularization prior instead (see prior_s_ba2()).
#
#   Rscript data-raw/estimate_s_ba2.R     # needs internet + nhanesA
# =============================================================================

suppressMessages(library(nhanesA))
for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) source(f)

trial_cols <- c("MGXH1T1","MGXH1T2","MGXH1T3","MGXH2T1","MGXH2T2","MGXH2T3")
fetch <- function(demo, body, grip) {
  d <- nhanes(demo)[, c("SEQN","RIAGENDR","RIDAGEYR")]
  b <- nhanes(body)[, c("SEQN","BMXWAIST","BMXHT")]
  g <- nhanes(grip)[, c("SEQN", trial_cols)]
  Reduce(function(x, y) merge(x, y, by = "SEQN"), list(d, b, g))
}
nh <- rbind(fetch("DEMO_G","BMX_G","MGX_G"), fetch("DEMO_H","BMX_H","MGX_H"))
nh$age <- nh$RIDAGEYR
nh$sex <- ifelse(nh$RIAGENDR %in% c("Male", 1), 1L, 2L)
tr     <- sapply(nh[trial_cols], as.numeric)
nh$grip <- apply(tr, 1, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
nh$whtr <- nh$BMXWAIST / nh$BMXHT
nh <- nh[nh$age >= 18 & nh$age <= 80 &
         complete.cases(nh[c("grip","whtr","age","sex")]), ]

coefs <- fit_kdm_coefficients(nh, c("grip","whtr"), source = "nhanes")

est_one <- function(sub, cf) {
  idx <- match(c("grip","whtr"), cf$bm)
  k <- cf$k[idx]; q <- cf$q[idx]; s <- cf$s[idx]
  den <- sum((k / s)^2)
  BAe <- ((sub$grip - q[1]) * (k[1]/s[1]^2) +
          (sub$whtr - q[2]) * (k[2]/s[2]^2)) / den
  list(s_ba2 = mean((BAe - sub$age)^2) - 1/den, estim_var = 1/den, n = nrow(sub))
}
for (sx in c(1, 2)) {
  r <- est_one(nh[nh$sex == sx, ], coefs[coefs$sex == sx, ])
  cat(sprintf("%s: n=%d  estimation-noise SD=%.0f yr  =>  s_BA^2=%.0f %s\n",
      ifelse(sx == 1, "Men  ", "Women"), r$n, sqrt(r$estim_var), r$s_ba2,
      ifelse(r$s_ba2 < 0, "(NEGATIVE -> degenerate, cannot estimate)", "")))
}
