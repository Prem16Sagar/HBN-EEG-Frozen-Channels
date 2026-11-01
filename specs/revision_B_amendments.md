# Amendments to revision_B_spec.md

Each amendment is recorded before the analysis it changes is run. Date: 22 October 2026.

## A1 — WP3: add a per-recording view (approved)

**Change.** In addition to the specified per-channel detection rates, report per recording: the number of sampled recordings with ≥ 1 frozen channel, how many had **all** frozen channels caught, how many had **≥ 1 missed** (Wilson 95% CIs), and a list of every recording with a miss (frozen count, missed count, bands of the missed channels).

**Reason.** Before the WP3 run finished, its progress log showed a single recording (R5 NDARWK983JFE) with 122 frozen channels of which 28 were caught. Per-channel rates would be dominated by that one recording. The per-recording view is added **alongside** the specified per-channel rates, never instead of them, and is labelled as added after this observation.

**Implementation.** A separate read-only cell reads WP3's saved `detector_xrel_per_channel.csv` and writes to a new folder. The WP3 runner is not modified.

## A2 — WP1: Releases 10 and 11 (approved conditionally)

**Facts established.** OpenNeuro lists Release 10 as ds005515 and Release 11 as ds005516, and its dataset descriptions give their FCP-INDI locations as `cmi_bids_R10` and `cmi_bids_R11` — the paths the census tried. The census reported only that `participants.tsv` was not found at those paths.

**Step 1 — diagnostic (read-only, no scanning).**
1. List the contents of `cmi_bids_R10/` and `cmi_bids_R11/` on FCP-INDI.
2. Check whether `participants.tsv` exists for ds005515 and ds005516 on OpenNeuro.
3. **Equivalence check:** download the resting-state `.set` file of three Release 3 children from both FCP-INDI and OpenNeuro (ds005507) and compare SHA-256 hashes. If Release 10 files exist on FCP-INDI, compare two of them across the sources as well.

**Step 2 — only if the diagnostic supports it.** Scan Releases 10 and 11 with the identical census rule, from whichever source the diagnostic validates, as a documented extension of census_v1. Any source other than FCP-INDI requires identical hashes in step 1.3; if the hashes differ, that difference is reported as a finding, and no data from the second source are mixed in. Step 2 will be specified in full, for approval, after the diagnostic.

## Also recorded: read-only checks of WP2 output

Two questions raised by WP2's results, answered by reading `harm_per_child.csv` only (no new analysis): (1) which child or children in the "whole-recording frozen only" group show the maximum inflation of 19.7×, with their frozen-channel counts; (2) which children lack an exponent under the wrong-order reference (n = 179 in the age-slope model, against 181 children).
