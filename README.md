# FitAge — Functional (Fitness) Age Calculator

An open-source R implementation that estimates a person's **functional age**
from physical performance measures, using the **Klemera-Doubal Method (KDM)** —
the same approach used for blood-biomarker biological age (e.g.
[`BioAge`](https://github.com/dayoonkwon/BioAge)), applied here to fitness markers.

Markers: grip strength, sit-to-stand / sit-and-rise, one-leg balance,
sit-and-reach flexibility, push-ups, waist-to-height ratio, reaction time.

## Why KDM (and not a single ML model)

A supervised model `age ~ f(grip, balance, …)` needs every marker measured on
the same people. No open dataset has all of these jointly. KDM removes that
requirement: **each marker is calibrated against age independently**, and the
markers are then combined. For marker *j* we fit

```
x_j = q_j + k_j · age      (residual SD = s_j)
```

and the functional-age estimate is

```
              Σ_j (x_j − q_j) · (k_j / s_j²)
FunctionalAge = ───────────────────────────────
                   Σ_j (k_j / s_j)²
```

The **weight of each marker falls out automatically** as `k_j / s_j²`: markers
that change steeply with age (`k` large) and are low-noise (`s` small) count
more. You never pick weights by hand, and markers can come from different
source datasets. With both `chrono_age` and a cohort-derived `s_ba2`, the
chronological-age-corrected KDM (`BA_EC`) is also supported.

## Data sources & status

NHANES 2011–2014 measures grip strength and waist+height on the **same**
individuals — that's the joint *anchor*, fitted directly. The remaining markers
are entered as published coefficients.

All non-NHANES markers are calibrated from **published summary tables** (means/SDs
or percentiles by age band), transcribed into tidy CSVs under `data-raw/sources/`
and fitted programmatically — see "How coefficients are fitted" below. UK Biobank
and CLSA individual data are access-controlled, so only their published norms are used.

| Marker | Source | Age range | Status |
|---|---|---|---|
| Grip strength (best single-hand) | NHANES 2011–2014 (max of 6 trials) | 18–80 | ✅ fitted (run build w/ internet) |
| Waist-to-height ratio | NHANES 2011–2014 (`BMXWAIST`/`BMXHT`) | 18–80 | ✅ fitted (run build w/ internet) |
| Push-ups (modified, 40 s) | Aandstad 2016, Norway (n=726) | 20–65 | ✅ fitted |
| Sit-and-reach flexibility | Aandstad 2016, Norway | 20–65 | ✅ fitted |
| One-leg balance | Aandstad 2016, Norway | 20–65 | ✅ fitted |
| Chair-rise / sit-to-stand | CLSA, Mayhew 2023 (Supp. App. 5) | 45–85 | ✅ fitted from exact published percentiles |
| Reaction time | UK Biobank (PMC8249619) | 40–70 | ⚠️ fitted, but residual SD is approximate |

> Notes: (1) Chair rise uses the exact per-year percentiles from Supplementary
> Appendix 5 (Tables 7/8), via `percentiles_to_moments()`.
> (2) **CLSA single-leg balance is deliberately not used** — it is censored at a
> 60 s ceiling and the authors could not even fit their GAMLSS model to it; the
> Norwegian balance test is used instead. (3) Reaction time's residual SD (`s`)
> is not reported by the source (~112 ms placeholder) and should be refined from
> the UK Biobank showcase; it carries little weight regardless. (4) Flexibility,
> balance and push-ups all come from the *same* Norwegian cohort, removing
> intercept mismatch between them. (5) A French max-push-up alternative
> (`data-raw/sources/pushups_french_nassif2012.csv`) is included but off by
> default — its female slope is unreliable. (6) Each marker is only valid within
> its calibration age range; applying it far outside that band extrapolates.

## How coefficients are fitted

Three ingestion paths, all producing the same `(k, q, s)` schema:

- **Joint microdata** (NHANES grip + waist-to-height) → `fit_kdm_coefficients()`
  regresses each marker on age per sex.
- **Grouped summary tables** (mean ± SD per age band) → `fit_from_grouped()`
  fits a weighted regression of band means on mid-age and pools within-band SDs.
- **Percentile tables** (e.g. CLSA P5…P95) → `percentiles_to_moments()` recovers
  per-band mean/SD assuming normality, then `fit_from_grouped()`.

## Honest limitations

- **Stitched populations.** Coefficients come from British, French, Canadian and
  US cohorts. Slope (`k`) and noise (`s`) transfer across populations reasonably;
  the level/intercept (`q`) does **not**, which biases an individual's score
  because `(x − q)` is taken against a different population's baseline. Fine for
  an educational/wellness tool — **not** a validated clinical biomarker.
- **Protocol matters.** Push-up cadence, sit-and-reach box offset, balance time
  cap, and reaction-time apparatus must match each source's protocol, or the
  calibration is invalid.
- **Conditional independence.** The simple KDM ignores cross-marker
  correlations (it uses the diagonal form, as `BioAge` does). A joint cohort
  would let you model the full covariance.
- **Weak markers need the chronological-age correction.** These fitness markers
  are individually weak age predictors (noise-to-slope of 24–109 years), so the
  *uncorrected* estimate (`BA_E`) is unstable with few markers and can return
  extreme or negative ages. See below.

## Corrected estimate (recommended)

The uncorrected Klemera-Doubal estimate has no anchor, so with few/weak markers
it swings wildly. The **corrected** estimate (`BA_EC`) adds a chronological-age
prior — a precision-weighted average of the marker estimate and your real age:

```r
functional_age(vals, female, chrono_age = 55, s_ba2 = prior_s_ba2(10))
```

`s_ba2` is the prior variance (how far functional age can sit from chronological
age). `prior_s_ba2(prior_sd_years = 10)` is the default — smaller `prior_sd`
pulls harder toward real age. A true KDM `s_ba2` needs a cohort with all markers
per person, which public data lacks for these measures, so this is a documented
modelling choice. With weak markers the prior carries real weight; with more/
stronger markers it automatically matters less. Example effect:

| Case | raw `BA_E` | corrected `BA_EC` |
|---|---|---|
| 55 yo woman, all 7 markers | 53.1 | 54.3 |
| 40 yo man, 3 markers | 17.6 | 35.4 |
| 40 yo man, 3 markers, very fit | −18.6 | 27.9 |

`try_fitage.R` uses the corrected estimate by default.

## Usage

```r
# load the package functions + bundled coefficients
for (f in list.files("R", full.names = TRUE)) source(f)
load("data/fitage_coefficients.rda")   # built by data-raw/build_coefficients.R

# a 55-year-old woman: pick the female calibration rows
female <- fitage_coefficients[fitage_coefficients$sex == 2, ]

vals <- c(grip = 28,            # best single-hand grip, kg
          whtr = 0.52,          # waist / height
          pushups_mod = 8,      # modified push-ups in 40 s
          balance_oneleg = 45,  # one-leg stand score (OLSsum)
          sit_reach = 20,       # sit-and-reach, cm
          chair_rise = 12.5,    # 5x sit-to-stand, seconds
          reaction = 540)       # reaction time, ms

res <- functional_age(vals, female, chrono_age = 55)

res$functional_age          #> 53.1   (estimated functional age)
res$functional_age_advance  #> -1.9   (functional minus chronological: ~2 yr "younger")
round(sort(res$contributions, decreasing = TRUE), 2)
#> balance_oneleg       reaction    pushups_mod     chair_rise           grip
#>          18.77           9.92           8.20           7.87           7.70
#>      sit_reach           whtr
#>           0.97          -0.37
```

You supply only the markers you have — KDM degrades gracefully with fewer.
By default, markers whose calibration does not cover the person's age are
dropped with a warning (e.g. applying the 45–85 chair-rise calibration to a
30-year-old); pass `enforce_age_band = FALSE` to override.

## Rebuilding the coefficient table

```bash
Rscript data-raw/build_coefficients.R   # NHANES step needs internet + the nhanesA package
```

## Tests

```bash
Rscript -e 'for (f in list.files("R", full.names=TRUE)) source(f); testthat::test_dir("tests/testthat")'
```

## License

MIT.
