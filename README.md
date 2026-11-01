# Frozen channels in the HBN-EEG releases

Code, specifications and derived results for:

> Prem Sagar and Ruban S. *Frozen channels are common in the Healthy Brain Network EEG releases and can bias average-referenced analyses.* Manuscript submitted to *Scientific Reports*.

The study scans every HBN-EEG resting-state recording in Releases 1–10 (2,597 recordings) for channels that hold one stored value for at least 1 s, measures what such channels do to an average reference and to the aperiodic exponent (Release 3, 181 children), tests whether PREP detects them, and audits how published HBN-EEG analyses order bad-channel detection and re-referencing.

**No EEG data are included.** All recordings are public: OpenNeuro ds005505–ds005515 (Releases 1–10) and the FCP-INDI mirror (`s3://fcp-indi/data/Projects/HBN/BIDS_EEG/`). The scripts download what they need.

## Repository layout

| Folder | Contents |
|---|---|
| `code/` | Every script, in one flat folder. Several scripts import others by name (`pilot_v31.py`, `wp1_census.py`), so keep them together and do not rename them. |
| `specs/` | Written specifications and amendments, each fixed before the analysis it governs was run. |
| `results/` | Derived tables reported in the paper (see below). |
| `figures/` | Figures 1-3 as PNG and PDF, plus `figure_numbers.json` with every number drawn. |

## How the scripts were run

All scripts were written as Google Colab cells with Google Drive mounted at `/content/drive`. Each one sets its paths at the top (`DRIVE = "/content/drive/MyDrive/hbn"`); change that line to run elsewhere, and put all of `code/` in that folder. Every run writes to a new time-stamped folder, refuses to overwrite earlier output, and stops if a reproduction gate fails. Exact package versions for each run are recorded in the provenance files the scripts write.

| Step | Script | Purpose | Paper |
|---|---|---|---|
| 0 | `setup_hbn_files.py` | Downloads and hash-checks three Release 3 pilot recordings | — |
| 1 | `pilot_v31.py` | Pipeline v3.1 (filtering, channel QC, average reference, Welch spectra, FOOOF 2–40 Hz). Importable module; run alone it processes the pilot | Methods |
| 2 | `run_r3_streaming.py` | Runs pipeline v3.1 on all Release 3 resting recordings, one child at a time | Methods |
| 3 | `analysis_r3_v31.py` | Pre-specified aperiodic analysis (plan fixed 20 Oct 2026; deviations D1, D2) | Suppl. Note 2 |
| 4 | `final_descriptives_cell.py` | Release 3 frozen-value types; attention coefficient intervals | Suppl. Table S1, Note 2 |
| 5 | `revision_P0_cell.py` | Decomposition of the fallback quantity, fit-range effects, exact stored values | Suppl. Table S2, Note 2 |
| 6 | `scan_stuck_channels.py` | Earlier 30-children-per-release scan; defines the 240-child cross-release sample | Gate G-C2; WP3 sample |
| 7 | `wp1_census.py` | Census of Releases 1–9 (gates G-C1, G-C2) | Results, Table 1 |
| 8 | `r10r11_diagnostic_cell.py` | Where Releases 10 and 11 can be read from (amendment A2, step 1) | Methods |
| 9 | `a2_r10r11_census.py` | Census of Releases 10 and 11 from OpenNeuro; combined summary (gates G-C3–G-C5) | Results, Table 1, Fig. 1 |
| 9b | `census_build_r10_tables.py` | Builds the Release 10 recording, channel and episode tables from the per-child results, using `wp1_census.py`'s own table code | Results tables |
| 10 | `wp2_harm.py` | Harm of frozen channels to the average reference, Release 3 (gate G-H1) | Results, Fig. 2, Suppl. Table S3 |
| 11 | `wp2_checks_cell.py` | Read-only checks on WP2 outputs (outliers, failed fits) | Results |
| 12 | `detector_test_runner.py` | PREP detection on Release 3 | Results, Suppl. Table S4 |
| 13 | `wp3_detector_xrel.py` | PREP detection on the cross-release sample, by episode duration (gate G-P1) | Results, Fig. 3, Suppl. Table S4 |
| 14 | `wp3_per_recording_cell.py` | Per-recording detection view (amendment A1) | Results, Suppl. Table S4 |
| 15 | `wp4_make_audit_tables.py` → `wp4_a3_update.py` → `wp4_spotcheck_update.py` | Pipeline audit tables: round 1, round 2 (amendment A3), author spot-checks. Run in that order in one folder; they reproduce `results/audit/` exactly | Results, Table 2, Suppl. Table S5 |
| 16 | `make_figures_v4.py` | Figures from the saved result files | Figs 1–3 |

## Results files

| File | What it holds |
|---|---|
| `results/census/combined_summary.csv` | Per-release and total counts with Wilson 95% CIs, Releases 1–11, with the selection rule for each row |
| `results/census/census_recordings_*.csv` | One row per scanned recording (Releases 1–9: 2,118; Release 10: 479): duration, frozen-channel counts, whether the 193,893.78125 value occurs |
| `results/census/census_channels_*.csv` | One row per frozen channel: location (scalp or peripheral), episodes, longest and total frozen time, duration band, stored values |
| `results/census/census_episodes_*.csv` | One row per frozen episode: channel, start, duration, stored value |
| `results/census/census_summary_R1-R9.csv`, `census_failures_R1-R9.csv` | Release 1–9 summary from `wp1_census.py`, and its scan failures (none) |
| `results/harm/harm_per_child.csv` | Per child: freeze group, exponent and reference amplitude under the correct, frozen-only and before-detection references, inflation and exponent error. Age and sex removed; they are in the HBN-EEG `participants.tsv` |
| `results/detection/detector_xrel_per_channel.csv` | Every frozen channel in the cross-release sample, its duration band and whether PREP flagged it |
| `results/audit/` | Screening log (97 records), audit table (12 studies) and the audit summary with its process notes |

Children with failed fits have empty exponent-error cells. Percentages in the paper use children with an estimable exponent as the denominator; `harm_summary.csv` from the original run used all children and is superseded on that point (Supplementary Table S3).

## Specifications

| File | Governs |
|---|---|
| `specs/revision_P0_spec.md` | `revision_P0_cell.py` |
| `specs/revision_B_spec.md` | WP1–WP4 |
| `specs/revision_B_amendments.md` | Amendments A1 and A2 |
| `specs/revision_B_A2_step2_spec_v2.md` | `a2_r10r11_census.py` |
| `specs/amendment_A3.md` | Second audit search round |

## Requirements

Python 3 with numpy, scipy, pandas, matplotlib, statsmodels, requests and mne, plus the two pinned packages in `requirements.txt` (`fooof==1.1.1`, `pyprep==0.9.0`).

## Licence and citation

Code is released under the MIT licence (`LICENSE`). The derived tables describe public HBN-EEG data; use of the underlying recordings is governed by their own terms. Please cite the paper and the archived release of this repository (`CITATION.cff`).
