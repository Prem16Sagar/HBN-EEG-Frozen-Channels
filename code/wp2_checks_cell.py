# =============================================================================
# WP2 READ-ONLY CHECKS (revision_B_amendments.md)      -- reads saved outputs only
#  1. Children in the "whole-recording frozen only" group with the largest inflation,
#     with their frozen scalp channels (from R3_v31_channels.csv).
#  2. Children without an exponent under the wrong-order reference (age-slope n = 179).
# Writes to a new folder revision_B/checks/wp2_checks_<timestamp>/.
# =============================================================================
import os, glob, time, numpy as np, pandas as pd
DRIVE = "/content/drive/MyDrive/hbn"
HARM = f"{DRIVE}/revision_B/harm/harm_v1"
OUT = f"{DRIVE}/revision_B/checks/wp2_checks_{time.strftime('%Y%m%d_%H%M%S')}"
summ = sorted(glob.glob(f"{HARM}/summary_*"))[-1]
H = pd.read_csv(f"{summ}/harm_per_child.csv")
C = pd.read_csv(f"{DRIVE}/R3_v31/R3_v31_channels.csv",
                usecols=["participant_id", "channel", "location", "invalid_reasons", "run_value", "longest_identical_run_s"])
P = pd.read_csv(f"{DRIVE}/R3_v31/R3_v31_participants.csv", usecols=["participant_id", "n_samples", "sfreq"])
os.makedirs(OUT, exist_ok=False)
W = H[H.group.str.startswith("whole")].sort_values("frozen_only__inflation_ratio", ascending=False)
cols = ["participant_id", "n_frozen_scalp_whole", "n_frozen_scalp_episodic", "n_frozen_peripheral",
        "frozen_only__inflation_ratio", "frozen_only__exponent_error", "wrong__inflation_ratio", "wrong__exponent_error"]
W[cols].to_csv(f"{OUT}/check1_whole_only_by_inflation.csv", index=False)
top = W[W.frozen_only__inflation_ratio > 1.1].participant_id.tolist()
F = C[C.participant_id.isin(top) & (C.location == "scalp") & C.invalid_reasons.fillna("").str.contains("stuck")]
F = F.merge(P, on="participant_id")
F["recording_s"] = F.n_samples / F.sfreq
F["samples_not_in_longest_run"] = np.round(F.n_samples - F.longest_identical_run_s * F.sfreq).astype(int)
F[["participant_id", "channel", "run_value", "longest_identical_run_s", "recording_s", "samples_not_in_longest_run",
   "invalid_reasons"]].to_csv(f"{OUT}/check1_frozen_channels_of_high_inflation_children.csv", index=False)
M = H[H[["correct__exponent", "wrong__exponent", "age", "sex"]].isna().any(axis=1)]
M[["participant_id", "group", "correct__exponent", "wrong__exponent", "correct__n_fit_fail", "wrong__n_fit_fail",
   "wrong__n_reference", "age", "sex"]].to_csv(f"{OUT}/check2_missing_from_age_slope.csv", index=False)
pd.set_option("display.width", 220)
print(f"source: {summ}/harm_per_child.csv\noutput: {OUT}\n")
print("CHECK 1 (check1_whole_only_by_inflation.csv): whole-recording-only children, top 5 by inflation")
print(pd.read_csv(f"{OUT}/check1_whole_only_by_inflation.csv").head(5).round(4).to_string(index=False))
print(f"\n  children in this group with frozen-only inflation > 1.1: {len(top)}")
print("\n  their frozen scalp channels (check1_frozen_channels_of_high_inflation_children.csv):")
print(pd.read_csv(f"{OUT}/check1_frozen_channels_of_high_inflation_children.csv").to_string(index=False))
print("\nCHECK 2 (check2_missing_from_age_slope.csv): children missing from the age-slope model")
print(pd.read_csv(f"{OUT}/check2_missing_from_age_slope.csv").round(4).to_string(index=False))
