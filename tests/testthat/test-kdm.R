# Synthetic round-trip: build markers with a known age relationship, fit
# coefficients, then confirm functional_age recovers age for an average person.

make_synthetic <- function(n = 4000, seed = 1) {
  set.seed(seed)
  age <- runif(n, 20, 80)
  sex <- rep(c(1, 2), length.out = n)
  data.frame(
    age  = age,
    sex  = sex,
    grip = 55 - 0.4 * age + rnorm(n, 0, 6),    # declines with age
    whtr = 0.40 + 0.0015 * age + rnorm(n, 0, 0.03)  # rises with age
  )
}

test_that("fit + functional_age recovers chronological age on average", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")

  # an "average 50-year-old male": values on the regression line
  male <- coefs[coefs$sex == 1, ]
  vals <- setNames((male$q + male$k * 50), male$bm)

  res <- functional_age(vals, male, chrono_age = 50)
  expect_equal(res$functional_age, 50, tolerance = 0.5)
  expect_equal(res$functional_age_advance, 0, tolerance = 0.5)
  expect_setequal(res$markers_used, c("grip", "whtr"))
})

test_that("stronger grip lowers functional age", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]

  on_line <- setNames(male$q + male$k * 50, male$bm)
  strong  <- on_line; strong["grip"] <- strong["grip"] + 10  # stronger than typical

  expect_lt(
    functional_age(strong, male)$functional_age,
    functional_age(on_line, male)$functional_age
  )
})

test_that("corrected estimate regularizes weak markers toward chronological age", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]

  # an extreme single-marker value that makes raw BA_E wild
  extreme <- c(grip = male$q[male$bm == "grip"] + male$k[male$bm == "grip"] * 50 + 40)

  raw  <- functional_age(extreme, male, chrono_age = 50)$functional_age
  corr <- functional_age(extreme, male, chrono_age = 50,
                         s_ba2 = prior_s_ba2(10))$functional_age

  # correction pulls the estimate back toward chronological age (50)
  expect_lt(abs(corr - 50), abs(raw - 50))
  # an unusually strong value still reads younger than 50, just not absurdly so
  expect_lt(corr, 50)
  expect_gt(corr, 0)
})

test_that("confidence band and stability reflect marker precision", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]
  vals <- setNames(male$q + male$k * 50, male$bm)

  res <- functional_age(vals, male, chrono_age = 50, s_ba2 = prior_s_ba2(10))
  expect_true(res$se > 0)
  expect_lt(res$ci95["low"], res$functional_age)
  expect_gt(res$ci95["high"], res$functional_age)
  expect_true(res$marker_information >= 0 && res$marker_information <= 1)
  # adding the prior tightens the band vs the uncorrected estimate
  raw <- functional_age(vals, male, chrono_age = 50)
  expect_lt(res$se, raw$se)
  # 2 markers -> "low" stability
  expect_equal(res$stability, "low")
})

test_that("leave_one_out ranks marker influence", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]
  vals <- c(grip = male$q[male$bm == "grip"] + male$k[male$bm == "grip"] * 50 + 8,
            whtr = male$q[male$bm == "whtr"] + male$k[male$bm == "whtr"] * 50)

  loo <- leave_one_out(vals, male, chrono_age = 50, s_ba2 = prior_s_ba2(10))
  expect_setequal(loo$marker, c("grip", "whtr"))
  # sorted by absolute influence; grip is the off-norm one so it moves more
  expect_gte(abs(loo$delta[1]), abs(loo$delta[2]))
})

test_that("prior_s_ba2 returns variance and validates input", {
  expect_equal(prior_s_ba2(10), 100)
  expect_error(prior_s_ba2(-1), "prior_sd_years > 0")
})

test_that("validate_coefficients rejects bad tables", {
  expect_error(validate_coefficients(data.frame(bm = "x")), "missing required")
  bad <- data.frame(bm = "x", sex = 1, k = 1, q = 0, s = -1)
  expect_error(validate_coefficients(bad), "must be positive")
})

test_that("age-band guard drops out-of-range markers", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]
  male$age_min <- 40; male$age_max <- 80   # calibration valid 40-80 only
  vals <- setNames(male$q + male$k * 50, male$bm)

  # age 30 is below the band -> both markers dropped -> error (nothing left)
  expect_error(suppressWarnings(functional_age(vals, male, chrono_age = 30)),
               "outside their calibration")
  # disabling the guard uses them anyway
  expect_silent(functional_age(vals, male, chrono_age = 30, enforce_age_band = FALSE))
  # in-band age is fine and reports nothing dropped
  res <- functional_age(vals, male, chrono_age = 60)
  expect_length(res$markers_dropped, 0)
})

test_that("missing markers degrade gracefully", {
  d <- make_synthetic()
  coefs <- fit_kdm_coefficients(d, c("grip", "whtr"), source = "synthetic")
  male <- coefs[coefs$sex == 1, ]
  res <- functional_age(c(grip = male$q[male$bm == "grip"] +
                            male$k[male$bm == "grip"] * 50), male)
  expect_equal(res$markers_used, "grip")
})
