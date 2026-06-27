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
