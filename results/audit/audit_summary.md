# WP4 — Pipeline audit: summary (run 26 October 2026)

## Procedure as run
- Five fixed queries, each run once, web search "standard" mode, top 10 results per query (50 records).
- Every record screened in search order; decision and reason logged (`audit_screening_log.csv`).
- Eligibility: original analysis of HBN EEG recordings with accessible methods. Excluded: descriptors, reviews, non-HBN data, non-research pages.
- Classification fields: data distribution, channel detection, detector type, reference, order category, location in paper (`audit_table.csv`).

## Yield
- 50 records -> 10 eligible (target was 30; the fixed queries were exhausted).
- 9 classified; 1 pending manual extraction (PMC12562131, full text blocked to Claude).

## Counts (n = 9 classified)
| Category | n |
|---|---|
| Detects before referencing | 3 (A04, A08, A09) |
| References before detecting | 2 (A02, A05) |
| No channel-level detection | 2 (A01, A06) |
| Not reported | 2 (A03, A07) |
| Pending | 1 (A10) |

## Descriptive observations (not judgements of any study's results)
1. Five different data distributions appear: original Biobank raw, HBN-provided preprocessed files, HBN-EEG BIDS (via EEGDash), the EEG2025 challenge preprocessed version (100 Hz), and a published post-ASR/ICA derivative. Four of nine studies used a pre-processed distribution, so upstream choices are inherited.
2. "Detects before referencing" does not guarantee a frozen channel is caught: A09 uses an upper RMS threshold only; A02's 5-s flat-line criterion cannot by itself flag shorter freezes. Detector type is recorded separately for this reason.
3. A02's reference is linked mastoids, not average, so its exposure depends only on the mastoid channels.

## Limitations
- n = 9 is small; counts are descriptive only. No percentages should be generalised to "the HBN literature".
- Search tool returns 10 results per query; query 4 (dataset identifiers) returned mostly noise.
- One exclusion (Q4 rank 1, medRxiv UM-EEG) was made from the record text without full-text check.

## Process error logged
While screening A05 (arXiv 2606.30104), Claude first quoted the paper's referencing sentence from memory with a citation index that pointed to a different document. The full text was then fetched; the sentence is present in Appendix A.1, and the classification now rests on that text. Logged here so the audit trail is complete.

## Spot-checks requested from the author (spec: author verifies a sample)
- A05: reading of the referencing sentence (Appendix A.1).
- A04: assignment to "Detects before referencing" when no voltage re-referencing occurs.
- A10: extract the four fields from PMC12562131.
- Q4 rank 1: confirm HBN is not used.

## Decision needed
The specification targeted 30 papers. Options:
- (a) Report n = 9-10 as is, with the limitation stated.
- (b) Amendment A3: re-run the same five queries in "extended" mode, screen new records by the same rules until 30 eligible or exhaustion. Decided after seeing these results, so it must be logged as a post hoc amendment motivated by sample size.

---
## Amendment A3 result (extended search, run 26 Oct 2026; see amendment_A3.md)
- 47 round-2 records screened (5 queries, extended mode). Mostly duplicates of round 1 and dataset/release pages.
- New eligible: 2. CoMET (arXiv 2509.00314) -> No channel-level detection. NeurML (ResearchGate) -> pending (blocked; one lookup search for an open copy found none).
- Two exclusions in round 2 made from record text without full-text check: arXiv 2011.12107, arXiv 2510.22364.
- Stopping rule reached: records exhausted. No further rounds (per A3 rule 5).

## Final counts (12 eligible; 10 classified, 2 pending)
| Category | n |
|---|---|
| Detects before referencing | 3 |
| No channel-level detection | 3 |
| References before detecting | 2 |
| Not reported | 2 |
| Pending manual extraction | 2 (A10 PMC12562131, A12 NeurML) |

Conclusion on sample: the fixed queries, in both modes, surface about a dozen eligible HBN EEG analyses. Report as counts with the search limitation stated; no percentages.

---
## Author spot-checks (26 Oct 2026) and final classification
- A05 (arXiv 2606.30104): author reads the sentence order as unclear -> reclassified Not determinable.
- A04 (arXiv 2503.19924): classification confirmed.
- A10 (PMC12562131), A12 (NeurML): not accessible to Claude or author -> reported as not accessible.
- Q4 rank 1 (medRxiv UM-EEG): author could not open it; exclusion rests on the record text only (stated as a limitation).

## FINAL (12 eligible)
| Category | n |
|---|---|
| Detects before referencing | 3 |
| No channel-level detection | 3 |
| References before detecting | 1 |
| Not reported | 2 |
| Not determinable | 1 |
| Not accessible | 2 |
