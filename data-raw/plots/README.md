# Plot-ready data for every FitAge marker

A single long-format CSV with the values needed to redraw the “marker vs age”
chart (points + fitted line + middle-50% band) for **every** marker in the
FitAge calibration table, plus a small slope-summary table for legends.

The numbers come from the same code path that builds
`data/fitage_coefficients.rda`, so the lines here are the lines FitAge
actually uses.

## Files

| File | What's in it |
| --- | --- |
| `plot_data_all.csv` | One long CSV: every marker × sex × (band points and fitted line). |
| `slope_all.csv` | One row per marker × sex: `k`, `q`, `s_pooled`, age range, source. |
| `export.R` | Regenerates both CSVs from `data-raw/sources/` + the bundled rda. |

Regenerate:

```sh
Rscript data-raw/plots/export.R
```

## Markers included

| `marker` | Label | Units | Source | Band points? |
| --- | --- | --- | --- | --- |
| `grip` | Grip strength | kg | NHANES 2011-2014 | – (individual data) |
| `whtr` | Waist-to-height ratio | ratio | NHANES 2011-2014 | – (individual data) |
| `sit_reach` | Sit-and-reach | cm | Kjær 2016 (Norway) | ✔ 5 bands × 2 sexes |
| `balance_oneleg` | One-leg balance (OLSsum) | s | Kjær 2016 (Norway) | ✔ 5 bands × 2 sexes |
| `pushups_mod` | Modified push-ups | reps / 40 s | Kjær 2016 (Norway) | ✔ 5 bands × 2 sexes |
| `sit_rise_floor` | Sit-rise floor test | score (0–10) | Araujo 2020 (provisional) | ✔ 8 bands × 2 sexes |
| `chair_rise` | 5-rep chair rise | s | Mayhew 2023 CLSA | ✔ 9 bands × 2 sexes (from percentiles) |
| `reaction` | Reaction time | ms | UK Biobank, derived | – (direct k/q/s only) |

French push-ups (Nassif 2012) is intentionally **off** in the FitAge model
(alternative protocol, female slope unreliable). Flip `include_french <- TRUE`
in `export.R` if you want to add it for comparison.

## `plot_data_all.csv` schema

Each row is **either** a published age-band data point (`row_type = "band"`)
**or** a point on the dense fitted line (`row_type = "fit"`). The marker-level
columns (`k`, `q`, `s_pooled`, `age_min`, `age_max`, `source`, …) repeat on
every row so a single filter on `marker` + `sex` gives you everything needed
to plot that subgroup.

| column | type | always present? | meaning |
| --- | --- | --- | --- |
| `marker` | id | ✔ | Marker key (e.g. `pushups_mod`) |
| `marker_label` | text | ✔ | Human-readable label |
| `units` | text | ✔ | Display units |
| `higher_is` | text | ✔ | `younger` or `older` (sign of the slope) |
| `sex_code` | int | ✔ | 1 = Men, 2 = Women |
| `sex` | text | ✔ | Same, labelled |
| `source` | text | ✔ | Provenance from the bundled coefficient table |
| `provisional` | bool | ✔ | `TRUE` for `sit_rise_floor` and `reaction` |
| `row_type` | enum | ✔ | `band` or `fit` |
| `age` | num | ✔ | Band midpoint (band rows) or grid age (fit rows) |
| `n_band` | int | band only | Participants in that age band |
| `mean` | num | band only | Group mean of the marker |
| `sd` | num | band only | Within-band SD |
| `p25`, `p75` | num | band only | mean ± 0.6745·SD (middle-50% endpoints) |
| `fit` | num | fit only | Line value: `q + k·age` |
| `lo50`, `hi50` | num | fit only | Ribbon: `fit ± 0.6745·SD(age)` |
| `k` | num | ✔ | Fitted slope (marker units per year) |
| `q` | num | ✔ | Fitted intercept |
| `s_pooled` | num | ✔ | Residual SD FitAge stores as `s` |
| `age_min`, `age_max` | num | ✔ | Calibration range (= line range) |
| `n_total` | int | ✔ | Sum of `n_band` (NA if no band data) |
| `n_bands` | int | ✔ | Number of age bands used for the fit (0 if none) |

### Ribbon SD

The `lo50`/`hi50` band uses a per-band SD interpolated linearly between the
published age midpoints (constant outside the band range). For markers
without per-band SDs (`grip`, `whtr`, `reaction`), the constant pooled
`s_pooled` is used, so the ribbon is parallel to the line.

## How to plot in any spreadsheet program

For one (`marker`, `sex`) pair, e.g. `pushups_mod` / Women:

1. **Filter** `plot_data_all.csv` to those rows.
2. **Points** — use rows with `row_type = "band"`: scatter `age` (x) vs `mean` (y);
   optionally use `p25`/`p75` as error-bar endpoints.
3. **Line** — use rows with `row_type = "fit"`: line plot of `age` vs `fit`.
4. **Ribbon** — use the same fit rows: filled band between `lo50` and `hi50`.
5. **Legend** — read `k`, `q`, `s_pooled`, `n_total`, `source` from any row
   (they're constant for the subgroup), or look them up in `slope_all.csv`.

To redraw the “modified push-ups by age” chart from the figure you supplied,
filter for `marker = "pushups_mod"` and plot Men and Women as separate series.
