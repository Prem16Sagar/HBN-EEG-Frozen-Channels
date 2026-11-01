# Amendment A2, step 2 — scanning Releases 10 and 11 (version 2, proposed, awaiting approval)

Date: 26 October 2026. Supersedes `revision_B_A2_step2_spec.md` (kept unchanged). Basis: `revision_B/checks/r10r11_diag_20261026_110501/diagnostic.json`.

## Changes from version 1 (self-review)

1. **Population from the folders that exist, not from participants.tsv.** Release 11's participants file lists 397 participants, while OpenNeuro's dataset descriptions give 430 subjects (Release 10: 533 in both). Selection now enumerates subject folders.
2. **Census code imported unchanged**, with only its download addresses redirected — no copied code.
3. **Wider equivalence gate (G-C3):** small sidecar files for ten children per release, plus one `.set` file per release.
4. **Version uncertainty recorded:** release descriptions and participants files hashed; the ambiguity between versions 1.0.0 and 1.0.1 for Release 10 stated, not resolved.
5. **Added descriptive split:** frozen-channel prevalence restricted to **scalp** channels (whole-recording vs episodic), for all releases, because harm via a v3.1-style reference concerns scalp channels only.

## Facts from the diagnostic

- No files exist at the documented FCP-INDI paths `cmi_bids_R10/` and `cmi_bids_R11/` (both list as empty). Other paths were not checked.
- OpenNeuro holds ds005515 (Release 10; participants file of 533 rows with RestingState, 480 "available" or "caution") and ds005516 (Release 11; 397 rows, no RestingState column).
- Three Release 3 resting-state `.set` files are byte-identical on FCP-INDI and OpenNeuro (ds005507), and match the hashes recorded for the pilot children.

## Source and provenance

OpenNeuro's public S3 copy (`https://s3.amazonaws.com/openneuro.org/<dataset>/…`). For ds005515 and ds005516, `dataset_description.json`, `participants.tsv` and (if present) `CHANGES` are saved with their SHA-256. The version served by S3 is not pinned; this is reported as a limitation.

## Gates (before any Release 10/11 scanning; failure stops the run)

- **G-C3 — source equivalence, Releases 1–9.** Per release: `channels.tsv` and `events.tsv` hashed on FCP-INDI and OpenNeuro for the first ten census participants (sorted), plus the resting-state `.set` of the first. All must be identical; any difference is reported with its file.
- **G-C4 — method identity.** The census module (`wp1_census.py`, imported unchanged from Drive) scans three Release 3 children with downloads redirected to OpenNeuro. Results must equal census_v1 exactly (frozen channels, episodes, values, longest run of every channel).

## Population

- **Primary rule (Releases 10 and 11):** every `sub-*` folder present on OpenNeuro is checked for a resting-state `.set`; children with the file are scanned. Children without it are recorded as "no resting-state file" and excluded from the denominator.
- **Release 10 cross-check against the Releases 1–9 rule:** report how many folders with a file are "available/caution" in its participants file, how many are listed otherwise or not listed, and how many "available/caution" participants lack a file. Release 10 prevalence is reported under both rules.
- **Release 11:** the participants file cannot be used for selection (no RestingState column); this is footnoted. Folders present but absent from the participants file are counted and reported.

## Method

The census rule, unchanged, through the imported census module: every run of identical consecutive samples lasting ≥ 1 s, on the full continuous recording, Cz exempt; the same duration bands.

## Outputs

- Per-child results in a new folder `census/census_v1_r10r11/per_child/`, written once and never rewritten.
- A **combined summary** in a new folder covering Releases 1–11, with the census tables (per release: the code, any frozen value, episodic frozen channels, Wilson CIs; dominant channels; episode durations; χ²), plus the **scalp-only split** of change 5 for every release. Releases 1–9 are read from census_v1 unchanged, and their gates G-C1 and G-C2 are carried forward. **Releases with different selection rules are reported separately and are not pooled into a single prevalence**; any all-release figure is given only as a census total with its rules stated.

## Runtime

About 1–1.5 hours (roughly 910 folders; OpenNeuro download speed untested), resumable.
