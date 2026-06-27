# FitAge — Methods & Data

How FitAge estimates a **functional ("fitness") age** from physical performance
measures, and where every number comes from.

---

## 1. The method: Klemera–Doubal (KDM)

FitAge **adapts the Klemera–Doubal framework to functional fitness markers.** KDM
was developed and validated for biological-age estimation from clinical/blood
biomarkers (e.g. the `BioAge` package); FitAge borrows the same mathematical
machinery but applies it to physical-performance measures. This is an
exploratory, educational use — fitness markers are individually *weaker* age
predictors than blood biomarkers, and FitAge has **not** been validated against
health outcomes the way blood-based clocks have. It is not equivalent to a
validated biological-age clock.

**Why KDM and not machine learning?** A supervised model `age ~ f(grip, balance,
…)` needs every marker measured on the *same individuals*. No open dataset has
all these fitness measures jointly. KDM removes that requirement: **each marker
is calibrated against chronological age independently**, then the markers are
combined. This is what lets us assemble markers from different source datasets.

**Per-marker calibration.** For each marker *j* (separately by sex) we fit a
linear age relationship and keep three numbers:

```
x_j = q_j + k_j · age      (residual SD = s_j)
```

- `k_j` — slope (how fast the marker changes per year of age)
- `q_j` — intercept
- `s_j` — residual standard deviation (scatter around the age line)

**Combination.** The marker-based functional age (`BA_E`) is:

```
              Σ_j (x_j − q_j) · (k_j / s_j²)
BA_E  =    ───────────────────────────────
                   Σ_j (k_j / s_j)²
```

Each marker's **weight is `k_j / s_j²`** — it counts more when it changes steeply
with age (large `k`) and is low-noise (small `s`). Weights are never set by hand.
The sign of `k` encodes direction automatically (grip declines, `k<0`; reaction
time rises, `k>0`).

**Chronological-age correction (recommended, default in the tooling).** These
fitness markers are individually weak age predictors, so `BA_E` is unstable with
few markers and can return extreme or negative ages. The corrected estimate
(`BA_EC`) adds a chronological-age prior — a precision-weighted average of the
marker estimate and the person's real age:

```
            Σ_j (x_j − q_j)(k_j/s_j²) + CA / s_ba2
BA_EC =   ──────────────────────────────────────────
                 Σ_j (k_j/s_j)² + 1 / s_ba2
```

**`s_ba2` here is a regularization choice, not a canonical KDM parameter
estimated from the assembled marker set.** In standard KDM, `s_ba2` is *estimated*
from a reference cohort with every marker measured per person. FitAge has no such
cohort. We did test estimating it from the one joint subset available — NHANES
grip + waist-to-height — and the KD estimator returns a **negative** variance
(men −96, women −141): those two markers are so weak that the estimation noise
(~31–33 yr SD) swamps the age signal, so the estimate is degenerate, not usable
(`data-raw/estimate_s_ba2.R`). So `s_ba2` is set deliberately as a prior
(`prior_s_ba2()`, default prior SD 10 yr) that regularizes the estimate toward
chronological age — not a quantity recovered from data. Smaller values pull
harder toward chronological age; tune to taste.

---

## 2. Turning published data into coefficients

Sources come in different shapes, so there are three ingestion paths — all
producing the same `(bm, sex, k, q, s)` schema (see `R/`):

| Path | Function | Used for |
|---|---|---|
| Joint individual microdata | `fit_kdm_coefficients()` | NHANES (grip, waist-to-height) — regress each marker on age per sex |
| Grouped means/SDs per age band | `fit_from_grouped()` | Norway, floor SRT — weighted fit of band means on mid-age; pool within-band SDs |
| Percentile tables | `percentiles_to_moments()` → `fit_from_grouped()` | CLSA chair rise — recover mean/SD per band assuming normality |
| Direct `k/q/s` | (used as-is) | Reaction time — slope/intercept derived analytically |

The whole table is rebuilt, fully auditably, by `data-raw/build_coefficients.R`
from the raw source files in `data-raw/sources/`.

---

## 3. The eight markers and their sources

| Marker (id) | Test / units | Source | Age range | How fitted |
|---|---|---|---|---|
| Grip strength (`grip`) | best single-hand, kg | **NHANES 2011–2014** (max of 6 per-hand trials) | 18–80 | joint microdata |
| Waist-to-height (`whtr`) | waist ÷ height | **NHANES 2011–2014** | 18–80 | joint microdata |
| Push-ups (`pushups_mod`) | modified (knee), reps in 40 s | **Aandstad 2016, Norway** (n=726) | 20–65 | grouped means/SDs |
| One-leg balance (`balance_oneleg`) | OLSsum score, max 120 s | **Aandstad 2016, Norway** | 20–65 | grouped means/SDs |
| Sit-and-reach (`sit_reach`) | cm | **Aandstad 2016, Norway** | 20–65 | grouped means/SDs |
| Chair rise (`chair_rise`) | 5× sit-to-stand, seconds | **CLSA / Mayhew 2023** (Supp. App. 5, exact percentiles) | 45–85 | percentile table |
| Floor Sitting-Rising (`sit_rise_floor`) | composite score 0–10 | **Araújo 2020** (n=6141) | 46–85 | grouped (median-based) ⚠️ provisional |
| Reaction time (`reaction`) | ms | **UK Biobank** (PMC8249619) | 40–70 | direct k/q/s ⚠️ SD approximate |

### Fitted coefficients (current build)

`x = q + k·age`, residual SD `s`; sex M=male, F=female.

| Marker | Sex | k | q | s | Age range |
|---|---|---|---|---|---|
| balance_oneleg | M | -0.8407 | 92.87 | 19.7 | 25–62.5 |
| balance_oneleg | F | -0.8686 | 96.51 | 21.2 | 25–62.5 |
| chair_rise | M | 0.0723 | 8.198 | 3.23 | 45–85 |
| chair_rise | F | 0.0887 | 7.37 | 3.30 | 45–85 |
| grip | M | -0.2355 | 56.32 | 8.79 | 18–80 |
| grip | F | -0.1572 | 36.41 | 5.68 | 18–80 |
| pushups_mod | M | -0.1602 | 19.13 | 4.12 | 25–62.5 |
| pushups_mod | F | -0.1255 | 13.66 | 4.04 | 25–62.5 |
| reaction | M | 2.85 | 361.6 | 112 | 40–70 |
| reaction | F | 3.40 | 345.2 | 112 | 40–70 |
| sit_reach | M | -0.1329 | 25.0 | 11.7 | 25–62.5 |
| sit_reach | F | -0.1022 | 28.87 | 13.3 | 25–62.5 |
| sit_rise_floor | M | -0.1159 | 13.67 | 2.8 | 48–83 |
| sit_rise_floor | F | -0.1427 | 15.20 | 2.8 | 48–83 |
| whtr | M | 0.00164 | 0.495 | 0.086 | 18–80 |
| whtr | F | 0.00145 | 0.534 | 0.103 | 18–80 |

### Marker strength (noise-to-slope, `s/|k|`, in years — lower = more weight)

`sit_rise_floor` 21.9 · `balance_oneleg` 23.9 · `pushups_mod` 28.9 ·
`reaction` 36.1 · `grip` 36.7 · `chair_rise` 40.9 · `whtr` 61.9 ·
`sit_reach` 109. Flexibility is the weakest age marker and contributes least —
exactly as KDM should weight it.

---

## 4. Key decisions

- **Sex-specific** calibration throughout (fitness norms differ strongly by sex).
- **Age-band guard** — `functional_age()` drops markers whose calibration does
  not cover the person's age, rather than extrapolating.
- **Grip = best single-hand** (max of the six NHANES per-hand trials), not the
  combined sum, to match home-dynamometer use and the CLSA dominant-hand norm.
- **Push-ups: Norwegian modified protocol** is the default; the French
  max-effort dataset is kept but disabled (its female slope is unreliable).
- **CLSA single-leg balance rejected** — censored at a 60 s ceiling; the authors
  could not fit their own GAMLSS model to it. Norwegian balance used instead.

## 5. Cross-checks performed

- **Synthetic round-trip** — an average person's values recover their
  chronological age exactly (unit test).
- **Grip** — NHANES single-hand fit (`q≈56`, `k≈−0.24`) lands near the CLSA
  dominant-hand chart (`k≈−0.29`).
- **Chair rise** — the exact-percentile fit matched an independent chart-read of
  the same paper's figure to two decimals.
- 15 automated tests (`tests/testthat/`).

## 6. Honest limitations

- **Stitched populations.** Markers come from US, Norwegian, Canadian, Brazilian
  and British cohorts. Slope (`k`) and noise (`s`) transfer reasonably across
  populations; the level/intercept (`q`) does **not**, which biases an
  individual's score. This is a wellness/educational tool — **not** a validated
  clinical biomarker.
- **Protocol dependence.** Each input must match its source's protocol (grip
  hand, push-up style, balance/reaction apparatus) or the calibration is invalid.
- **Conditional independence.** KDM ignores cross-marker correlation (diagonal
  form). `chair_rise` and `sit_rise_floor` overlap (both "getting up"), so that
  capacity is mildly over-weighted; we cannot quantify it without a joint cohort.
- **Provisional pieces.** The floor SRT is median-based with an approximate SD
  (Brazilian clinical sample); the reaction-time SD (~112 ms) is approximate.
- **Weak markers / age coverage.** Several markers cover only older adults, and
  all are individually weak age predictors — hence the chronological-age
  correction is recommended.

---

## Sources

- **NHANES 2011–2014** Muscle Strength & Body Measures — CDC/NCHS (public).
- Aandstad A, et al. *Normative values for musculoskeletal- and neuromotor
  fitness in apparently healthy Norwegian adults.* BMC Sports Sci Med Rehabil
  2016 (PMC5116214).
- Mayhew AJ, et al. *Normative values for grip strength, gait speed, TUG, single
  leg balance, and chair rise from the CLSA.* Age and Ageing 2023; 52(4):afad054.
- Araújo CGS, et al. *Sitting–rising test: sex- and age-reference scores derived
  from 6141 adults.* Eur J Prev Cardiol 2020; 27(8):888.
- Reaction time vs age in UK Biobank — PMC8249619.
