# Analysis specification — revision items P0-1, P0-2, P0-4

**Status:** proposed, awaiting approval. No code has been written or run.
**Date:** 22 October 2026
**Scope:** Revision_Roadmap.md items P0-1 (decompose D), P0-2 (fit-range sensitivity) and P0-4 (exact stored frozen values). Nothing else.

---

## 0. Rules

1. **Read-only on existing material.** No existing code, data, result or manuscript file is modified or overwritten.
2. **New output folder.** All outputs go to a new, timestamped folder: `/content/drive/MyDrive/hbn/revision_P0/run_<YYYYmmdd_HHMMSS>/`. The script refuses to run if that folder already exists.
3. **Every reported number is read from a file this run produces**, and the file is named when the number is reported.
4. **No new model fits.** All quantities are recomputed from the per-channel fits already stored by pipeline v3.1. No new fit ranges, R² thresholds, PHF definitions or criterion-2 variants.
5. **Estimation only.** No new hypothesis tests or p-values. Uncertainty is reported as bootstrap confidence intervals.
6. **Not run:** anything under "Do not do" in Revision_Roadmap.md — no Q2 or Q3, no new PHF or criterion-2 variants, no multiverse of fit settings or censoring thresholds, no additional detectors.

## 1. Inputs (all existing, on Drive)

| File | Used for |
|---|---|
| `R3_v31/analysis/analysis_sample.csv` | Defines the analysis sample (181 children) and supplies saved participant-level summaries, age, sex, attention, frontal HF |
| `R3_v31/R3_v31_channels.csv` | Per-channel fits (`ap2_exponent`, `ap2_r2`, `ap1_exponent`, `ap1_r2`), channel status (`location`, `is_invalid`, `is_high_variance`, `in_reference`), frozen values (`invalid_reasons`, `run_value`, `longest_identical_run_s`) |
| `R3_v31/R3_v31_participants.csv` | `sfreq`, `n_samples` (to classify whole-recording vs episodic) |
| `stuck_scan/per_participant/R{1,2,4,5,6,7,8,9}_*.csv` | Frozen values in the cross-release sample. Files prefixed `R3_` (the scanner's validation re-scan) are **excluded** to avoid double counting Release 3 |

The SHA-256 of every input file is recorded in the run log.

## 2. Common conventions

- **Sample:** the 181 children in `analysis_sample.csv`. The attention analyses use complete cases (expected n = 177).
- **Bootstrap:** 2,000 resamples of children with replacement; percentile 95% CIs; `numpy.random.default_rng(20261022)` — a new seed, distinct from the original analysis (20261020). Differences between two quantities use **paired** resampling (same child indices for both).
- **Regression:** Huber M-estimator (statsmodels `RLM`, `HuberT`), as in the pre-specified analysis. Formulas are evaluated in an isolated namespace. **Coefficients are selected by name**, never by position.
- **Failed fits:** a channel whose fit raised an error has NaN exponent and NaN R². As in pipeline v3.1, NaN exponents are excluded from every summary, and a NaN R² counts as failing the R² filter. The number of NaN fits per child is reported.
- **Units:** exponent differences in raw exponent units, and also divided by the between-child SD of the standard exponent (`scalp__ap2_exponent`), stated as "SD units".

---

## 3. P0-1 — Decompose D

### 3.1 Channel sets and summaries (per child, 2–40 Hz fit)

- **V** (all valid scalp) = channels with `location == "scalp"` and not `is_invalid`.
- **R** (reference-eligible) = V and `in_reference` (equivalently, V and not `is_high_variance`).
- **F** (fit filter) = `ap2_r2 ≥ 0.85`.

Eight summaries form a 2 × 2 × 2 factorial — statistic {median, mean} × channel set {R, V} × filter {F, none}. Notation `stat(set, filter)`.

### 3.2 Gate G1 — reproduction check (run must stop if it fails)

Four of the eight summaries correspond to saved columns and must reproduce them for every child (absolute difference ≤ 1 × 10⁻⁹):

| Recomputed | Saved column |
|---|---|
| S = median(R, F) | `scalp__ap2_exponent` |
| B = median(R, none) | `scalp__ap2_exponent_nofilter` |
| A = median(V, none) | `scalp_allvalid__ap2_exponent_nofilter` |
| U = mean(V, none) | `scalp_allvalid__ap2_exponent_nofilter_mean` |

D = U − S must also reproduce the saved `D` column. **If any check fails, the run stops and reports the mismatch; no decomposition is produced.**

### 3.3 Primary decomposition (order fixed in advance — the reviewer's order)

| Term | Definition | Isolates |
|---|---|---|
| C | B − S | removing the R² filter (censoring) |
| H | A − B | re-including high-variance channels |
| G | U − A | mean instead of median (aggregation) |

**Gate G2:** C + H + G = D for every child (absolute difference ≤ 1 × 10⁻¹²).

For D and each term, report: median across children with bootstrap 95% CI; mean; SD; and the percentage of children for whom the term is exactly zero (|term| ≤ 1 × 10⁻¹²).

### 3.4 Order sensitivity

Sequential decompositions depend on order. All 3! = 6 orderings of the three factors (filter, channel set, statistic) are computed from the eight factorial cells. For each factor, report its contribution under each ordering (median across children, bootstrap CI), and its **Shapley contribution** — the per-child average over the six orderings — which is order-independent. The primary order (3.3) is reported first; the others are sensitivity results.

### 3.5 Subgroups

- **Nothing excluded:** zero high-variance channels in V **and** zero channels in V failing F (NaN R² counts as failing).
- **At least one excluded:** all other children.

Report n and median D with bootstrap CI for each group. **Gate G3:** in the "nothing excluded" group, D must equal G exactly (|D − G| ≤ 1 × 10⁻¹²), confirming that D reduces to mean − median there.

---

## 4. P0-2 — Fit-range sensitivity

### 4.1 Exponents compared

- **Matched set (primary):** for each child, the median over R channels passing R² ≥ 0.85 in **both** fits (`ap2_r2 ≥ 0.85` and `ap1_r2 ≥ 0.85`), for each fit range: E2m (2–40 Hz) and E1m (1–40 Hz). This isolates the fit range, because Task 1 found that the as-reported shift also differs in channel set.
- **As reported (secondary):** E2 = `scalp__ap2_exponent`, E1 = `scalp__ap1_exponent` — the quantities behind Figure 3b, each with its own R² filter.

The shift is Δ = E1 − E2, computed for both versions. The number of children whose matched set is empty is reported; they are excluded from the matched-set analyses.

### 4.2 Statistics

1. **Distribution of Δ:** median, mean, SD, IQR, 5th and 95th percentiles; also in SD units.
2. **Agreement between fit ranges:** Spearman ρ(E1, E2); ICC(3,1) (two-way, consistency, single measures, as requested by the reviewer) and ICC(2,1) (two-way, absolute agreement), both from the Shrout–Fleiss two-way ANOVA formulas; bootstrap CIs for all three.
3. **Does Δ track predictors?** Spearman ρ between Δ and age, log₁₀ frontal relative 70–110 Hz power, and attention; bootstrap CIs.
4. **Do effects change with fit range?** Each model is fitted under both fit ranges, and the difference is estimated with a paired bootstrap CI:
   - **Age slope:** exponent (raw units) ~ age (years, linear) + C(sex), Huber. Slope in exponent units per year; also standardised (× SD of age ÷ SD of exponent).
   - **Attention coefficient:** exactly the pre-specified fallback model — z(exponent) ~ z(attention) + cr(age, df = 3, centred) + C(sex), Huber, coefficient selected by name; n = 177.

### 4.3 Decision rule for the manuscript claim (fixed now, before any result)

The claim that the fit-range choice **changes conclusions** is retained only if, under the matched-set primary analysis, the paired 95% CI for the change in the age slope **or** the change in the attention coefficient excludes zero. Otherwise the manuscript uses the narrower wording: *"aperiodic summaries are sensitive to fit range and summary statistic; these must be reported."* In either case, the median-ratio comparison ("17 times") is dropped, and shifts are reported in SD units.

---

## 5. P0-4 — Exact stored frozen values

1. **Frozen channels.** Release 3: rows of `R3_v31_channels.csv` whose `invalid_reasons` contains "stuck". Other releases: rows with `stuck == True` in the scan files.
2. **Whole vs episodic.** Release 3: `longest_identical_run_s ≥ n_samples / sfreq − 1/sfreq` (as in the descriptives cell). Scan: the saved `whole_recording` column.
3. **For every distinct `run_value`,** report: the value at full precision (`repr`, 17 significant digits); its float32 hexadecimal representation and whether it is exactly representable in float32; whether its exact negative also occurs (mirror check); and counts of channels, children, releases, whole vs episodic, and scalp vs peripheral (peripheral = the 19 channels in pipeline v3.1).
4. **Specific answer to Inventory §5.1:** state whether values of magnitude ≈ 220,853 occur in Release 3, with the exact value or values.
5. **No sample extremes are computed.** The draft's "≈ ±199,700" statement was a comparison of different quantities (Task 1) and will be removed rather than re-measured.

---

## 6. Outputs (new run folder)

| File | Contents |
|---|---|
| `run_log.json` | Input paths and SHA-256, package versions (python, numpy, scipy, pandas, statsmodels, patsy), seed, n_boot, gate results G1–G3 |
| `P0-1_per_child.csv` | Eight factorial summaries, D, C, H, G, the Shapley contributions, NaN-fit and exclusion counts, subgroup label |
| `P0-1_summary.csv` | Medians, CIs, means, SDs and %-zero for D and each term, all six orderings and Shapley |
| `P0-1_subgroups.csv` | D by subgroup |
| `P0-2_per_child.csv` | E1, E2, E1m, E2m, Δ (both versions), matched-set channel count |
| `P0-2_results.json` | Distribution, agreement, Δ–predictor correlations, slopes and coefficients under each fit range with paired-difference CIs, and the decision-rule outcome |
| `P0-4_stuck_values.csv` | One row per distinct frozen value with all fields in §5.3 |

## 7. Execution and reporting

- One Colab cell, run by you with Drive mounted. Expected runtime: a few minutes (the bootstraps fit Huber models only; no EEG processing).
- **Before you run it,** I will test the cell on synthetic files with known answers, including deliberately failing each gate, and in a polluted notebook namespace.
- **After the run,** you send me the printed summary. I report each number with the output file it came from. If a result contradicts the review, I report it plainly.
