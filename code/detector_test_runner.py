# =============================================================================
# DETECTOR TEST: do standard bad-channel tools catch the stuck channels?   (Colab)
#
# Runs PREP's noisy-channel detection (pyprep NoisyChannels, default criteria) on
# every processed Release 3 recording, and scores it against the stuck channels
# identified by pipeline v3.1 (ground truth: identical-value run >= 1 s).
#
# Design (fixed before running):
#   * Detection on the FULL continuous recording, stored values (uV -> V), as a
#     typical user would run it; pyprep's own detrending (1 Hz) applied.
#   * Default criteria WITHOUT RANSAC for every child. RANSAC can only ADD flags, so
#     it is run ONLY for children where a stuck channel was missed (adaptive).
#   * Scored per stuck channel: caught or not, by which criterion; split by
#     whole-recording vs episodic, and by stuck-value type. Cz (reference) ignored.
#   * Context: how many NON-stuck scalp channels each child also had flagged.
#   * Resumable + crash-safe: one result file per child, written atomically.
# Run after: drive.mount('/content/drive')
# =============================================================================
import os, sys, json, time, shutil, tempfile, subprocess, gc, warnings
import numpy as np, pandas as pd, requests, scipy.io

try:
    import mne, pyprep
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "mne", "pyprep==0.9.0"])
    import mne, pyprep
from pyprep import NoisyChannels
mne.set_log_level("ERROR"); warnings.filterwarnings("ignore")

DRIVE = "/content/drive/MyDrive/hbn"
R3_DIR = f"{DRIVE}/R3_v31"
OUT = f"{DRIVE}/detector_test"
S3 = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG/cmi_bids_R3"
SEED = 20261020
MAX_NEW = None                      # children to process this run (None = all remaining)
KNOWN, SAT = 193893.78125, 220853.4
PERIPH = {"E48","E119","E127","E126","E113","E49","E56","E107","E99","E63",
          "E94","E68","E43","E120","E88","E73","E81","E125","E128"}


def atomic_json(obj, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(path + ".tmp", path)


def fetch(url, dest, retries=3):
    for k in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(dest + ".part", "wb") as f:
                    for b in r.iter_content(1 << 20): f.write(b)
            os.replace(dest + ".part", dest); return
        except Exception as e:
            if k == retries: raise
            time.sleep(2 * k)


def ground_truth():
    C = pd.read_csv(os.path.join(R3_DIR, "R3_v31_channels.csv"),
                    usecols=["participant_id", "channel", "location", "invalid_reasons",
                             "run_value", "longest_identical_run_s"])
    P = pd.read_csv(os.path.join(R3_DIR, "R3_v31_participants.csv"), usecols=["participant_id", "sfreq", "n_samples"])
    C = C.merge(P, on="participant_id", how="left")
    C["stuck"] = C.invalid_reasons.fillna("").str.contains("stuck")
    C["whole"] = C.longest_identical_run_s >= C.n_samples / C.sfreq - 1.0 / C.sfreq
    C["type"] = np.where(np.isclose(C.run_value, KNOWN, atol=.05), "193,893.78 code",
                np.where(np.isclose(C.run_value.abs(), SAT, atol=1.0), "±220,853", "other constant"))
    return C


def detect(raw, ransac):
    nc = NoisyChannels(raw, do_detrend=True, random_state=SEED)
    nc.find_all_bads(ransac=ransac)
    return {k: sorted(str(c) for c in v) for k, v in nc.get_bads(as_dict=True).items()}


def run_child(pid, truth):
    base = f"{S3}/sub-{pid}/eeg/sub-{pid}_task-RestingState"
    tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
    try:
        sp, cp = os.path.join(tmp, "eeg.set"), os.path.join(tmp, "channels.tsv")
        fetch(base + "_eeg.set", sp); fetch(base + "_channels.tsv", cp)
        m = scipy.io.loadmat(sp, squeeze_me=True, struct_as_record=False)
        names = pd.read_csv(cp, sep="\t")["name"].tolist()
        raw = mne.io.RawArray(np.asarray(m["data"], np.float64) * 1e-6,
                              mne.create_info(names, float(m["srate"]), "eeg"))
        raw.set_montage(mne.channels.make_standard_montage("GSN-HydroCel-129"), on_missing="ignore")
        t0 = time.time(); bads = detect(raw, ransac=False); used_ransac = False
        stuck = truth[truth.stuck]
        missed = [c for c in stuck.channel if c not in bads["bad_all"]]
        if missed:                                   # adaptive: RANSAC only adds flags
            bads = detect(raw, ransac=True); used_ransac = True
        rows = []
        for _, r in stuck.iterrows():
            crit = [k for k, v in bads.items() if k != "bad_all" and r.channel in v]
            rows.append(dict(channel=r.channel, location=r.location, type=r.type, whole=bool(r.whole),
                             caught=r.channel in bads["bad_all"], criteria=";".join(crit),
                             caught_without_ransac=bool(r.channel not in missed)))
        scalp_ok = [c for c in truth[(truth.location == "scalp") & ~truth.stuck].channel]
        return dict(participant_id=pid, seconds=round(time.time() - t0, 1), used_ransac=used_ransac,
                    n_stuck=len(stuck), stuck=rows,
                    n_nonstuck_scalp=len(scalp_ok),
                    n_nonstuck_scalp_flagged=int(sum(c in bads["bad_all"] for c in scalp_ok)),
                    flags_by_criterion={k: len(v) for k, v in bads.items()})
    finally:
        shutil.rmtree(tmp, ignore_errors=True); gc.collect()


def summarise():
    d = os.path.join(OUT, "per_participant")
    J = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".json") and not f.endswith(".failed.json"):
            with open(os.path.join(d, f)) as fh: J.append(json.load(fh))
    S = pd.DataFrame([dict(participant_id=j["participant_id"], **s) for j in J for s in j["stuck"]])
    print(f"\nDETECTOR TEST: {len(J)} children, {len(S)} stuck channels")
    if len(S):
        def rate(g): return f"{g.caught.sum()} of {len(g)} ({100 * g.caught.mean():.1f}%)"
        print(f"  caught overall                 {rate(S)}")
        print(f"  caught without RANSAC          {S.caught_without_ransac.sum()} of {len(S)} ({100 * S.caught_without_ransac.mean():.1f}%)")
        for lab, g in (("whole recording", S[S.whole]), ("episodic", S[~S.whole])):
            if len(g): print(f"  {lab:<30} {rate(g)}")
        for t, g in S.groupby("type"): print(f"  type {t:<25} {rate(g)}")
        for loc, g in S.groupby("location"): print(f"  location {loc:<21} {rate(g)}")
        crit = S.criteria.str.split(";").explode().replace("", np.nan).dropna().value_counts()
        print("  criteria that caught them: " + ", ".join(f"{k} ({v})" for k, v in crit.items()))
        miss = S[~S.caught]
        if len(miss):
            print("  MISSED: " + ", ".join(f"{r.participant_id}:{r.channel} ({'whole' if r.whole else 'episodic'}, {r.type})"
                                         for _, r in miss.head(15).iterrows()))
        S.to_csv(os.path.join(OUT, "detector_stuck_channels.csv"), index=False)
    nf = pd.Series([j["n_nonstuck_scalp_flagged"] for j in J])
    print(f"  context: non-stuck scalp channels also flagged per child: median {nf.median():.0f}, "
          f"IQR {nf.quantile(.25):.0f}-{nf.quantile(.75):.0f}, max {nf.max():.0f}")
    print(f"  RANSAC needed for {sum(j['used_ransac'] for j in J)} children")
    fails = [f for f in os.listdir(d) if f.endswith(".failed.json")]
    print(f"  failed children: {len(fails)}")


def main():
    if OUT.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    d = os.path.join(OUT, "per_participant"); os.makedirs(d, exist_ok=True)
    G = ground_truth()
    pids = sorted(G.participant_id.unique())
    todo = [p for p in pids if not os.path.exists(os.path.join(d, f"{p}.json"))
            and not os.path.exists(os.path.join(d, f"{p}.failed.json"))]
    print(f"pyprep {pyprep.__version__}, mne {mne.__version__} | {len(pids)} children, "
          f"{len(pids) - len(todo)} done, {len(todo)} to process")
    if MAX_NEW is not None: todo = todo[:MAX_NEW]
    for k, pid in enumerate(todo, 1):
        t0 = time.time()
        try:
            res = run_child(pid, G[G.participant_id == pid])
            atomic_json(res, os.path.join(d, f"{pid}.json"))
            caught = sum(s["caught"] for s in res["stuck"])
            msg = f"OK  stuck {res['n_stuck']}, caught {caught}{' (RANSAC used)' if res['used_ransac'] else ''}"
        except Exception as e:
            atomic_json({"participant_id": pid, "reason": f"{type(e).__name__}: {e}"}, os.path.join(d, f"{pid}.failed.json"))
            msg = f"FAIL {type(e).__name__}: {e}"
        print(f"[{k}/{len(todo)}] {pid}  {msg}  ({time.time() - t0:.0f} s)")
    summarise()


if __name__ == "__main__":
    main()
