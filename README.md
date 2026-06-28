# FitAge — Functional (Fitness) Age Calculator

An open-source R implementation that estimates a person's **functional age**
from physical performance measures. FitAge adapts the Klemera-Doubal Method
(KDM) framework to functional fitness markers.
It is an **exploratory, educational** tool: fitness markers are
individually weaker age predictors than blood biomarkers, and FitAge has **not**
been validated against health outcomes — it is **not** equivalent to a validated
biological-age clock.

Markers: 
- Grip strength
- Sit-to-stand / sit-and-rise
- One-leg balance
- Sit-and-reach flexibility
- Push-ups
- Waist-to-height ratio
- Reaction time
- Floor Sitting-Rising Test

Created and maintained by **Filipa Santos Rodrigues**.

> **Related:** [TrueAge](https://trueage.thehumanruntime.com) is the
> blood-biomarker sibling of FitAge — the same Klemera-Doubal backbone, but a
> different input layer (blood chemistry instead of physical-performance tests).

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
| Grip strength (best single-hand) | NHANES 2011–2014 (max of 6 trials) | 18–80 | ✅ fitted |
| Waist-to-height ratio | NHANES 2011–2014 (`BMXWAIST`/`BMXHT`) | 18–80 | ✅ fitted |
| Push-ups (modified, 40 s) | Kjær 2016, Norway (n=726) | 20–65 | ✅ fitted |
| Sit-and-reach flexibility | Kjær 2016, Norway | 20–65 | ✅ fitted |
| One-leg balance | Kjær 2016, Norway | 20–65 | ✅ fitted |
| Chair-rise / sit-to-stand | CLSA, Mayhew 2023 (Supp. App. 5) | 45–85 | ✅ fitted from exact published percentiles |
| Floor Sitting-Rising Test (0–10) | Araújo 2020 (n=6141) | 46–85 | ⚠️ provisional — median-based, SD approximate |
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
> (7) The floor Sitting-Rising Test is **provisional**: calibrated on the median
> composite score per age band (the score is discrete and ceiling-censored below
> ~45, so it is restricted to 46+), with an IQR-derived pooled SD. Brazilian
> clinical sample (Araújo 2020) — refine with exact P25/P75 from the source.

## How coefficients are fitted

Three ingestion paths, all producing the same `(k, q, s)` schema:

- **Joint microdata** (NHANES grip + waist-to-height) → `fit_kdm_coefficients()`
  regresses each marker on age per sex.
- **Grouped summary tables** (mean ± SD per age band) → `fit_from_grouped()`
  fits a weighted regression of band means on mid-age and pools within-band SDs.
- **Percentile tables** (e.g. CLSA P5…P95) → `percentiles_to_moments()` recovers
  per-band mean/SD assuming normality, then `fit_from_grouped()`.

## Honest limitations

- **Stitched populations.** Coefficients come from US (NHANES), Norwegian (Kjær),
  Canadian (CLSA), Brazilian (Araújo) and British (UK Biobank) cohorts.
  Slope (`k`) and noise (`s`) transfer across populations reasonably;
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

## Corrected estimate 

The uncorrected Klemera-Doubal estimate has no anchor, so with few/weak markers
it swings wildly. The **corrected** estimate (`BA_EC`) adds a chronological-age
prior — a precision-weighted average of the marker estimate and your real age:

```r
functional_age(vals, female, chrono_age = 55, s_ba2 = prior_s_ba2(10))
```

`s_ba2` here is a **regularization choice, not a canonical KDM parameter
estimated from the marker set.** Standard KDM *estimates* `s_ba2` from a cohort
with every marker measured per person; we have none, and the one joint subset
available (NHANES grip + waist-to-height) makes the KD estimator degenerate — it
returns a negative variance because those markers are too weak (see
`data-raw/estimate_s_ba2.R`). So `s_ba2` is set as a deliberate prior:
`prior_s_ba2(prior_sd_years = 10)` is the default, and smaller `prior_sd` pulls
harder toward real age. With weak markers the prior carries real weight; with
more/stronger markers it automatically matters less. Example effect:

| Case | raw `BA_E` | corrected `BA_EC` |
|---|---|---|
| 55 yo woman, all 8 markers | 54.5 | 54.8 |
| 40 yo man, 3 markers | 17.6 | 35.4 |
| 40 yo man, 3 markers, very fit | −18.6 | 27.9 |

`try_fitage.R` uses the corrected estimate by default.

## Usage

```r
# load the package functions + bundled coefficients
for (f in list.files("R", full.names = TRUE)) source(f)
load("data/fitage_coefficients.rda")   # built by data-raw/build_coefficients.R

# a 35-year-old woman: pick the female calibration rows
female <- fitage_coefficients[fitage_coefficients$sex == 2, ]

# at age 35, chair rise (45-85), floor sit-rise (48-85) and reaction time
# (40-70) aren't calibrated, so supply the markers that apply at this age:
vals <- c(grip = 26,            # best single-hand grip, kg
          pushups_mod = 20,     # modified push-ups in 40 s
          balance_oneleg = 60,  # one-leg stand score (OLSsum)
          sit_reach = 10,       # sit-and-reach, cm
          whtr = 0.41)          # waist / height

# corrected estimate (recommended)
res <- functional_age(vals, female, chrono_age = 35, s_ba2 = prior_s_ba2(10))

res$functional_age      #> 30.5   (estimated functional age)
res$ci95                #> low 13.7, high 47.2   (95% band; SE = res$se = 8.6)
res$marker_information  #> 0.27   (27% of the estimate is from markers, 73% the age prior)
res$stability           #> "moderate"  (5 markers, 0 provisional)

# how much does each marker move the score?
leave_one_out(vals, female, chrono_age = 35, s_ba2 = prior_s_ba2(10))
#>           marker provisional fa_without  delta
#> 1    pushups_mod       FALSE       36.6  +6.16
#> 2           grip       FALSE       28.3  -2.12
#> 3           whtr       FALSE       32.2  +1.71
#> 4 balance_oneleg       FALSE       28.8  -1.62
#> 5      sit_reach       FALSE       29.8  -0.67
```

You supply only the markers you have — KDM degrades gracefully with fewer.
Every score comes with a **confidence band**, a **stability** flag, and
`marker_information` (how much of the estimate is real markers vs the age prior),
plus `leave_one_out()` to see each marker's influence. By default, markers whose
calibration does not cover the person's age are dropped with a warning (e.g. the
45–85 chair-rise calibration on a 30-year-old); pass `enforce_age_band = FALSE`
to override.

## Rebuilding the coefficient table

```bash
Rscript data-raw/build_coefficients.R   
```

## Tests

```bash
Rscript -e 'for (f in list.files("R", full.names=TRUE)) source(f); testthat::test_dir("tests/testthat")'
```

## Author & citation

Created by **Filipa Santos Rodrigues**. The methodology (assembling the
coefficient table from public sources and adapting the Klemera-Doubal framework
to functional fitness markers) is original work by the author.

If you use FitAge, please cite it — see [`CITATION.cff`](CITATION.cff), or:

> Santos Rodrigues, F. (2026). *FitAge: Functional (Fitness) Age Calculator* (v0.1.0). https://github.com/FilipaCSR/FitAge

## License

MIT © 2026 Filipa Santos Rodrigues. See [`LICENSE`](LICENSE).
