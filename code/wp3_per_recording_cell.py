# =============================================================================
# WP3 PER-RECORDING VIEW (amendment A1)                -- reads saved outputs only
# Added after one recording (122 frozen channels, 28 caught) appeared in the WP3 log;
# reported ALONGSIDE the specified per-channel rates, never instead of them.
# Writes to a new folder revision_B/detector_xrel/per_recording_<timestamp>/.
# =============================================================================
import os, glob, time, numpy as np, pandas as pd
DRIVE = "/content/drive/MyDrive/hbn"
XREL = f"{DRIVE}/revision_B/detector_xrel/xrel_v1"
EARLIER = f"{DRIVE}/detector_test/detector_stuck_channels.csv"
OUT = f"{DRIVE}/revision_B/detector_xrel/per_recording_{time.strftime('%Y%m%d_%H%M%S')}"

def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return f"[{100 * (c - h):.1f}, {100 * (c + h):.1f}]"

summ = [s for s in sorted(glob.glob(f"{XREL}/summary_*")) if os.path.exists(f"{s}/detector_xrel_per_channel.csv")][-1]
D = pd.read_csv(f"{summ}/detector_xrel_per_channel.csv")
os.makedirs(OUT, exist_ok=False)
def per_rec(frame, label):
    g = frame.groupby(["participant_id"] + (["release"] if "release" in frame else []))
    R = g.agg(n_frozen=("caught", "size"), n_caught=("caught", "sum")).reset_index()
    R["n_missed"] = R.n_frozen - R.n_caught
    if "band" in frame:
        mb = frame[~frame.caught].groupby("participant_id").band.agg(lambda b: "; ".join(f"{k}: {v}" for k, v in b.value_counts().items()))
        R["missed_bands"] = R.participant_id.map(mb).fillna("")
    R.insert(0, "source", label); return R
R = per_rec(D, "cross-release sample")
R.to_csv(f"{OUT}/per_recording.csv", index=False)
rows = []
for label, frame in (("cross-release sample", R),) + ((("Release 3 (earlier run)", per_rec(pd.read_csv(EARLIER), "Release 3 (earlier run)")),) if os.path.exists(EARLIER) else ()):
    n = len(frame); m = int((frame.n_missed > 0).sum())
    rows.append(dict(source=label, recordings_with_frozen=n, all_caught=n - m, with_missed=m,
                     pct_with_missed=100 * m / n, ci_with_missed=wilson(m, n),
                     missed_channels=int(frame.n_missed.sum()), frozen_channels=int(frame.n_frozen.sum())))
pd.DataFrame(rows).to_csv(f"{OUT}/per_recording_summary.csv", index=False)
pd.set_option("display.width", 220)
print(f"source: {summ}/detector_xrel_per_channel.csv\noutput: {OUT}\n")
print("PER-RECORDING SUMMARY (per_recording_summary.csv):")
print(pd.read_csv(f"{OUT}/per_recording_summary.csv").round(1).to_string(index=False))
miss = pd.read_csv(f"{OUT}/per_recording.csv"); miss = miss[miss.n_missed > 0]
print(f"\nRECORDINGS WITH >= 1 MISSED CHANNEL (per_recording.csv): {len(miss)}")
if len(miss): print(miss.to_string(index=False))
