# WP2b specification: controls and realistic references for the harm analysis
# (amendment A4 to revision_B_spec.md)

**Status:** proposed, awaiting approval. No code exists yet.
**Written:** 1 November 2026, after the second external review (items 17, 19–21, 24–26, 28) and after the exploratory analyses in `review2_exploratory.py` had been seen. WP2b is therefore a post hoc follow-up, and will be reported as one. What this document fixes in advance is how WP2b is computed and how each result will be interpreted, before any WP2b output exists.

## 1. Question

WP2 showed that adding partly frozen scalp channels to an average reference inflates scalp amplitude and shifts the aperiodic exponent. WP2b asks:

1. **Transients or dilution?** Is the harm caused by the frozen channels' content (step transients after filtering), or only by the extra channels in the average (dilution)?
2. **Do realistic references protect?** Is a median reference, computed without any detection, robust to frozen channels? Does interpolating flagged channels before averaging behave like excluding them?
3. **Is the harm measure-dependent?** Inflation over 2–40 Hz (the band of the exponent fit) as well as 1–110 Hz; exponent error on a fixed channel set; alpha power as a second measure.
4. **Does timing matter?** Freeze episodes and transitions inside the analysed eyes-closed segments, rather than over the whole recording.
5. **What does "fit failure" mean?** Separate fitting errors from fits below R² 0.85.

## 2. Data and processing (unchanged from WP2)

- The 181 Release 3 children analysed by pipeline v3.1; QC flags read unchanged from `R3_v31_channels.csv`.
- Processing imported unchanged from `pilot_v31.py`: `load_set`, `ec_blocks`, `preprocess_channel` (1 Hz high-pass and notches per channel, before referencing), Welch per eyes-closed block (2 s Hann, 50% overlap), `fit_aperiodic` (FOOOF 1.1.1, fixed, 2–40 Hz), `spectral_features`.
- Summary set for the exponent, valid-scalp set for amplitude, and the R² ≥ 0.85 rule: as in WP2.

## 3. Reference conditions

| Condition | Channels averaged | Status |
|---|---|---|
| `correct` | valid, non-high-variance scalp + Cz | as WP2 (gate G-H1) |
| `wrong` | all scalp + Cz (before detection) | as WP2 (gate G-H2) |
| `frozen_only` | correct + frozen scalp channels | as WP2 (gate G-H2) |
| `zero_k` | correct + k constant-zero channels, where k = number of frozen scalp channels added in `frozen_only` | **new: dilution control** |
| `median` | channel-wise median over all scalp + Cz, no detection | **new: robust reference** |
| `interp_avg` | v3.1-invalid scalp channels replaced by spherical-spline interpolation (MNE, GSN-HydroCel-129 positions), then the average of all scalp + Cz | **new: interpolate-then-average** |

`zero_k` equals `correct` scaled by N/(N+k), so it changes only the denominator; the difference between `frozen_only` and `zero_k` is the effect of the frozen channels' content. A reference built from the channels PREP's detection leaves in is not needed: in Release 3, PREP flagged all 348 frozen channels, so such a reference contains no frozen channel by construction. This is stated as a reason, not tested.

## 4. Measures, per child and condition

1. Inflation: median valid-scalp RMS divided by its value under `correct`, over 1–110 Hz (as WP2) and over 2–40 Hz (new; same integration and line-bin exclusion).
2. Exponent: median over summary channels with R² ≥ 0.85 (as WP2), and the error relative to `correct`.
3. Fixed-set exponent error: the median over channels that reach R² ≥ 0.85 under both `correct` and the condition (new).
4. Alpha: median over valid scalp channels of 8–13 Hz absolute power, as a ratio to `correct` (new).
5. Fit outcomes per condition: channels with a fitting error (NaN), and channels fitted with R² < 0.85, counted separately (new).

## 5. Freeze timing (per child, new)

From the raw stored data, using the census rule (identical value for ≥ 1 s, Cz exempt): for each frozen scalp channel, the number of episodes, transitions (freeze onsets plus offsets) and frozen seconds, over the whole recording and inside the eyes-closed segments.

## 6. Gates (any failure stops the run)

- **G-H1:** the `correct` exponent reproduces the saved v3.1 exponent within 10⁻⁶ for every child (as WP2).
- **G-H2:** the `wrong` and `frozen_only` exponents and 1–110 Hz inflation reproduce WP2's `harm_per_child.csv` within 10⁻⁶ for every child.
- **G-H3:** full-recording episode counts per frozen scalp channel equal the census (`census_episodes_R1-R9.csv`, Release 3) for every channel.
- Synthetic tests before any real data: `zero_k` with k = 0 equals `correct` exactly; a planted step-freeze inflates `frozen_only` but not `zero_k`; interpolation of an intact channel reproduces it closely; transition counts match planted episodes.

## 7. Pre-stated analyses and how each result will be read

Groups as in WP2 (partly frozen, n = 55; frozen throughout only, n = 45; none, n = 81). Medians with exact order-statistic 95% CIs; percentages of children with |error| > 0.05 use children with an estimable exponent as the denominator, with fitting failures reported separately.

| # | Comparison | Reading fixed now |
|---|---|---|
| 1 | `zero_k` vs `frozen_only`, partly frozen group | If median `zero_k` inflation is below 1.1-fold while `frozen_only` is above 2-fold, the harm is attributed to content (transients), not dilution. If `zero_k` inflation is at least half of `frozen_only` on the log scale, dilution is reported as a substantial contributor. Anything between is reported as mixed. |
| 2 | `zero_k` vs `frozen_only` for NDARJP444DZM (28 channels) | Reported as a single case, whichever way it falls; no generalisation. |
| 3 | `median` vs `correct` in each group | `median` is called protective if, in the partly frozen group, median inflation is below 1.1-fold and fewer than 10% of children have |error| > 0.05. Its cost in children without frozen channels is reported alongside, whatever it is. |
| 4 | `interp_avg` vs `correct` in each group | Reported as the agreement between two ways of excluding the same channels; no threshold. |
| 5 | 2–40 Hz vs 1–110 Hz inflation | Both reported; the headline uses whichever band the paper's claim concerns, with the other alongside. |
| 6 | Fixed-set vs WP2 exponent error | Both reported; if the medians differ by more than 0.05, the text says that part of the WP2 error reflects a change of channel set. |
| 7 | Alpha ratio | Reported as a second measure; no threshold. |
| 8 | Within-segment dose–response, partly frozen group | Spearman ρ of log inflation and |error| with transitions inside the eyes-closed segments, frozen seconds inside them, and number of partly frozen channels. Reported as exploratory. |
| 9 | The 18 children with inflation below 1.1-fold | How many have no freeze transition inside the eyes-closed segments; reported as counts. |

Whatever the outcomes, they replace the corresponding exploratory statements in the manuscript, and the Discussion's mechanism paragraph is rewritten to match them.

## 8. Outputs

New folder `revision_B/harm/harm_v2b/`, resumable per child, never overwriting: `per_child/*.json`; `summary_<timestamp>/` with `harm2b_per_child.csv`, `harm2b_timing_per_channel.csv`, `harm2b_summary.csv` (every comparison in §7), `gates_provenance.json` (gates, package versions, seed). The printed summary names the file behind every number.

## 9. Execution

- One Colab cell, `wp2b_harm.py`, importing `pilot_v31.py` unchanged; `mne` needed for interpolation.
- Expected runtime: similar to WP2 (each recording is downloaded once; the extra conditions add Welch and FOOOF fits only), roughly 1–2 hours.
- Before you run it, I test it on synthetic data with known answers (§6) and in a fresh namespace. After the run, you send the printed summary; I report each number with its source file.

## 10. Not included (and why)

- A second release or a connectivity measure (review item 28): reported as a limitation; alpha power is the added measure.
- Relaxed R² refits for the two unfittable children (item 16): the worst-case imputation already bounds their effect.
- PREP with line-noise removal (item 33) and MNE's flat annotation as a comparator (item 32): optional, separate from WP2b.
