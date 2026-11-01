# Analysis specification — Option B work package (P1-1, P1-2, P1-3, P1-4)

**Status:** proposed, awaiting approval. No code written or run.
**Date:** 22 October 2026
**Decisions in force:** D-1 = B (data-integrity-led); D-2 = full census; D-3 = descriptive pipeline audit.

---

## 0. Rules (unchanged from the P0 work)

1. Read-only on all existing code, data, results and the manuscript.
2. Each work package writes only to its own new, timestamped folder under `/content/drive/MyDrive/hbn/revision_B/`; runs refuse to overwrite.
3. Every reported number is read from a file the run produces, and the file is named.
4. Estimation only: bootstrap CIs, no new p-values (except the χ² homogeneity statistic, already part of the draft).
5. Nothing on the roadmap's "Do not do" list: no Q2/Q3, no new PHF or criterion-2 variants, no multiverse of fit settings, at most one detector (PREP, as already used).
6. Every work package is tested on synthetic data with known answers, including its gates, before you run it.

## Order and runtime

| WP | Item | What | Runtime (Colab) |
|---|---|---|---|
| WP1 | P1-3 | Census of frozen channels in every resting-state recording, all releases, with episode durations | ~3–4 h, resumable |
| WP2 | P1-1 | Cohort-level harm in Release 3: reference inflation and exponent error under wrong-order referencing | ~3–4 h, resumable |
| WP3 | P1-4 | PREP detection on the cross-release sample, by episode duration | ~1.5 h, resumable |
| WP4 | P1-2 | Structured audit of published HBN-EEG preprocessing | Desk work (done by me with web search; you spot-check) |

WP1 runs first: WP2 and WP3 use its episode data.

---

## WP1 — Census (P1-3)

**Population.** Every participant with RestingState status "available" or "caution" in each release's `participants.tsv`, at `cmi_bids_R{1..11}` on the FCP-INDI S3 bucket. Releases R10 and R11 are attempted; if a release's `participants.tsv` is not found, it is reported as absent, not assumed. Release 3 is **re-scanned** with the census code so that all releases use one identical method.

**Per channel (full continuous recording, stored values, no filtering):**
- every maximal run of identical consecutive samples lasting ≥ 1 s (the v3.1 rule), with each run's start, duration and exact value;
- derived: frozen (≥ 1 such run; Cz exempt), number of episodes, total frozen duration, longest run, whole-recording (single run spanning every sample), distinct frozen values;
- non-finite sample count.

**Duration bands (fixed now).** Each frozen channel is assigned by its **longest** run: [1, 5) s; [5, 30) s; [30 s, whole recording); whole recording.

**Gates.**
- **G-C1:** for Release 3, census frozen verdicts must equal pipeline v3.1's verdicts on every channel (`R3_v31_channels.csv`), and its longest-run durations must match within one sample.
- **G-C2:** for the 240 children in the earlier cross-release sample, verdicts must equal the earlier scan files.
Any mismatch stops the summary step and is reported.

**Outputs.**
- `census_channels.csv` — one row per frozen channel: release, participant, channel, location (scalp/peripheral from the v3.1 channel list), episodes, total and longest frozen duration, whole-recording flag, duration band, frozen value(s) at full precision. **This is the per-recording list proposed for release with the paper.**
- `census_recordings.csv` — one row per recording: number of frozen channels by type and band; any-episodic flag.
- `census_summary.csv` — per release: n recordings scanned and failed, n and % with ≥ 1 channel at 193,893.78125, at any frozen value, and with ≥ 1 **episodic** frozen channel, each with Wilson 95% CIs; dominant channels; episode-duration distribution (median, IQR, by band).
- An all-releases total **computed over the census itself**, so it needs no weighting. χ² homogeneity across releases reported as before.

## WP2 — Cohort-level harm in Release 3 (P1-1)

**Population.** The 181 analysed Release 3 children.

**Processing.** Identical to pipeline v3.1 (1 Hz high-pass, line-noise notches, eyes-closed blocks, Welch, fooof 1.1.1 fixed 2–40 Hz) — imported from the existing `pilot_v31.py`, not copied or edited. The **only** thing varied is which channels form the average reference. v3.1's quality-control flags are used unchanged for everything else, including which channels enter the exponent summaries.

| Reference | Channels averaged | Purpose |
|---|---|---|
| **Correct** (v3.1) | valid, non-high-variance scalp + Cz | baseline |
| **Wrong order** | **all** scalp channels + Cz (frozen and high-variance included) | what a pipeline does if it references before detecting |
| **Frozen-only** | correct set + frozen scalp channels | isolates the contribution of frozen channels alone |

**Gate G-H1:** the "correct" condition must reproduce each child's saved `scalp__ap2_exponent` within 1 × 10⁻⁶. Mismatch stops the run.

**Per child, for each of the two contaminated references:**
- **Inflation ratio** = median over valid scalp channels of 1–110 Hz RMS (line-harmonic bins excluded), contaminated reference ÷ correct reference.
- **Exponent error** = contaminated `scalp__ap2_exponent` − correct.
- **Fit failures**: number of valid scalp channels failing R² ≥ 0.85.

**Groups.** Children with no frozen scalp channel; whole-recording frozen scalp channels only; ≥ 1 episodic frozen scalp channel. (Frozen *peripheral* channels never enter the reference and are counted separately.)

**Statistics.** Per group and reference: median and IQR of the inflation ratio and exponent error, bootstrap CI of the median (2,000 resamples, seed 20261023), and the percentage of children with |exponent error| > 0.05. **Consequence for an association (secondary, fixed now):** age slope under the wrong-order reference versus the correct one, same Huber model as P0-2 (exponent ~ age + C(sex)), paired bootstrap CI of the difference.

**Outputs.** `harm_per_child.csv`, `harm_summary.csv`, `harm_results.json`.

## WP3 — PREP detection across releases (P1-4)

**Population.** The 240 children of the earlier cross-release sample (fixed sample files reused), plus the existing Release 3 detector results.

**Detector.** Exactly as in the Release 3 test: pyprep 0.9.0 `NoisyChannels`, default criteria, `do_detrend=True`, full continuous recording in volts; RANSAC run only for a child whose frozen channel was missed.

**Ground truth.** WP1 census rows for those children.

**Reported.** Detection rate overall, by duration band, by frozen-value type and by release, each with Wilson 95% CIs; the criteria that caught each channel; the number of missed channels with their bands and values. Release 3 results are shown alongside, labelled as from the earlier run. The caveat that line-noise removal was not applied is carried forward.

**Output.** `detector_xrel_per_channel.csv`, `detector_xrel_summary.csv`.

## WP4 — Audit of published HBN-EEG preprocessing (P1-2)

**Search (fixed queries, run once, logged with date).** Web search for: "Healthy Brain Network" EEG preprocessing; "HBN-EEG"; "HBN" EEG resting state; the OpenNeuro identifiers ds005505–ds005516; "EEG Foundation Challenge" HBN. Sources: PubMed/PMC, bioRxiv/medRxiv, arXiv, OpenReview, journal pages.

**Eligibility.** Original analyses of HBN EEG recordings (any task), 2017–2026, with accessible full text or methods. Excluded: reviews, dataset descriptors, and papers whose methods cannot be accessed.

**Selection rule (to prevent cherry-picking).** Records are screened in the order the fixed queries return them; the **first 30 eligible** are included. Every screened record and its include/exclude reason is logged.

**Extraction per paper.** Data distribution used (original HBN Biobank files / HBN-EEG BIDS releases / distributed preprocessed version); whether and how bad channels are detected; the order of detection relative to re-referencing; reference used. Recorded as a short paraphrase with section location — no long quotations.

**Classification (fixed categories).** Detects before referencing · References before detecting · No channel-level detection · Not reported · Not determinable.

**Verification.** Each classification is checked a second time against the source text; you spot-check a random third.

**Reporting.** Counts per category only, plus the table in the supplement. **No statement about whether any paper's results were affected.**

**Output.** `audit_screening_log.csv`, `audit_table.csv`, `audit_summary.md`.

**Deferred.** Checking whether the distributed preprocessed (challenge) version still contains frozen channels needs those files located first, and needs a different detection rule, because filtering turns a constant into a near-constant rather than identical samples. It will be specified separately once the files are found.

---

## Not in this package

Curator contact (D-4) — your earlier decision stands; the reviewer recommends it, and it remains open. The origin checks (P2-1), simulation robustness (P2-3) and all writing changes come after this package.
