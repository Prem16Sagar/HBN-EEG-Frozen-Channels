# =============================================================================
# WP3 -- PREP DETECTION ACROSS RELEASES, BY EPISODE DURATION     (Colab cell)
# Implements revision_B_spec.md, WP3 (approved).
#
# * Population: the 240 children of the earlier cross-release sample (sample_R*.csv).
# * Detector: pyprep 0.9.0 NoisyChannels, default criteria, do_detrend=True, full continuous
#   recording in volts; RANSAC only for a child whose frozen channel was missed. Identical to
#   the Release 3 test. No line-noise removal (caveat carried forward).
# * Ground truth: frozen channels from the WP1 census (latest summary whose gates passed).
# * Gate G-P1: detection re-run on 3 Release 3 children from the earlier detector test must
#   reproduce their caught/criteria results exactly.
# Resumable: per-child results written once; summaries to new folders.
# Run after: drive.mount('/content/drive'). WP1 census must be complete.
# =============================================================================
import os, sys, json, time, glob, shutil, tempfile, subprocess, gc, warnings
import numpy as np, pandas as pd, requests, scipy.io

try:
    import mne, pyprep
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "mne", "pyprep==0.9.0"])
    import mne, pyprep
from pyprep import NoisyChannels
mne.set_log_level("ERROR"); warnings.filterwarnings("ignore")

DRIVE   = "/content/drive/MyDrive/hbn"
CENSUS  = f"{DRIVE}/revision_B/census/census_v1"
SAMPLE  = f"{DRIVE}/stuck_scan"                          # sample_R*.csv from the earlier scan
EARLIER = f"{DRIVE}/detector_test"                       # Release 3 detector test (per_participant/*.json)
RUN     = f"{DRIVE}/revision_B/detector_xrel/xrel_v1"
S3_ROOT = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG"
SEED, MAX_NEW, G_P1_N = 20261020, None, 3
KNOWN = 193893.78125
BANDS = ["[1, 5) s", "[5, 30) s", "[30 s, whole)", "whole recording"]


class GateFailure(RuntimeError):
    pass


def fetch(url, dest, retries=3):
    last = None
    for k in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                if r.status_code == 404: raise FileNotFoundError(f"404: {url}")
                r.raise_for_status()
                with open(dest + ".part", "wb") as f:
                    for b in r.iter_content(1 << 20): f.write(b)
            os.replace(dest + ".part", dest); return
        except FileNotFoundError:
            raise
        except Exception as e:
            last = e; time.sleep(2 * k)
    raise IOError(f"failed after {retries} attempts: {last}")


def load_any_set(set_path, base_url):
    m = scipy.io.loadmat(set_path, squeeze_me=True, struct_as_record=False)
    if "data" in m:
        data, sf, nb, pn = m["data"], float(m["srate"]), m.get("nbchan"), m.get("pnts")
    elif "EEG" in m:
        E = m["EEG"]; data, sf, nb, pn = E.data, float(E.srate), E.nbchan, E.pnts
    else:
        raise RuntimeError("unsupported .set layout")
    if isinstance(data, str):
        fdt = os.path.join(os.path.dirname(set_path), data)
        fetch(base_url.rsplit("/", 1)[0] + "/" + data, fdt)
        data = np.fromfile(fdt, dtype=np.float32).reshape(int(pn), int(nb)).T
    return np.asarray(data, dtype=np.float64), sf


def detect(raw, ransac):
    nc = NoisyChannels(raw, do_detrend=True, random_state=SEED)
    nc.find_all_bads(ransac=ransac)
    return {k: sorted(str(c) for c in v) for k, v in nc.get_bads(as_dict=True).items()}


def run_child(rel, pid, frozen):
    base = f"{S3_ROOT}/cmi_bids_{rel}/sub-{pid}/eeg/sub-{pid}_task-RestingState"
    tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
    try:
        sp, cp = os.path.join(tmp, "eeg.set"), os.path.join(tmp, "channels.tsv")
        fetch(base + "_eeg.set", sp); fetch(base + "_channels.tsv", cp)
        X, sf = load_any_set(sp, base + "_eeg.set")
        names = pd.read_csv(cp, sep="\t")["name"].tolist()
        raw = mne.io.RawArray(X * 1e-6, mne.create_info(names, sf, "eeg"))
        raw.set_montage(mne.channels.make_standard_montage("GSN-HydroCel-129"), on_missing="ignore")
        bads = detect(raw, ransac=False); used_ransac = False
        missed = [c for c in frozen.channel if c not in bads["bad_all"]]
        if missed:
            bads = detect(raw, ransac=True); used_ransac = True
        rows = [dict(channel=r.channel, caught=r.channel in bads["bad_all"],
                     caught_without_ransac=r.channel not in missed,
                     criteria=";".join(k for k, v in bads.items() if k != "bad_all" and r.channel in v))
                for _, r in frozen.iterrows()]
        return dict(release=rel, participant_id=pid, used_ransac=used_ransac, rows=rows, bads=bads)
    finally:
        shutil.rmtree(tmp, ignore_errors=True); gc.collect()


def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def vtype(values):
    vals = [float(v) for v in str(values).split(";") if v]
    if KNOWN in vals: return "193,893.78 code"
    if any(abs(abs(v) - 220853.4) < 1 for v in vals): return "±220,853"
    if -193894.0 in vals: return "-193,894"
    return "other constant"


def census_channels():
    ok = []
    for s in sorted(glob.glob(f"{CENSUS}/summary_*")):
        g, c = os.path.join(s, "gates.json"), os.path.join(s, "census_channels.csv")
        if os.path.exists(g) and os.path.exists(c):
            with open(g) as fh: gates = json.load(fh)
            if all(v.get("passed") is True for v in gates.values()):
                ok.append(s)
    if not ok:
        raise RuntimeError("no census summary with all gates passed -- run WP1 to completion first")
    return ok[-1], pd.read_csv(os.path.join(ok[-1], "census_channels.csv"))


def main():
    if RUN.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    d = os.path.join(RUN, "per_child"); os.makedirs(d, exist_ok=True)
    csum, CC = census_channels()
    print(f"pyprep {pyprep.__version__}, mne {mne.__version__} | ground truth: {csum}")

    # ---- Gate G-P1: reproduce the earlier Release 3 detector results exactly
    prior = sorted(glob.glob(f"{EARLIER}/per_participant/*.json"))
    prior = [p for p in prior if not p.endswith(".failed.json")]
    chosen = []
    for p in prior:
        with open(p) as fh: j = json.load(fh)
        if j["n_stuck"] > 0: chosen.append(j)
        if len(chosen) == G_P1_N: break
    g1 = []
    for j in chosen:
        fz = CC[(CC.release == "R3") & (CC.participant_id == j["participant_id"])]
        res = run_child("R3", j["participant_id"], fz)
        old = {s["channel"]: (s["caught"], s["criteria"]) for s in j["stuck"]}
        new = {r["channel"]: (r["caught"], r["criteria"]) for r in res["rows"]}
        g1.append(dict(participant_id=j["participant_id"], identical=old == new, old=old, new=new))
    gate = dict(children=len(g1), all_identical=all(x["identical"] for x in g1), detail=g1)
    print(f"GATE G-P1: {len(g1)} Release 3 children re-run; identical to earlier results: {gate['all_identical']}")
    if not gate["all_identical"] or not g1:
        S = os.path.join(RUN, f"summary_{time.strftime('%Y%m%d_%H%M%S')}"); os.makedirs(S, exist_ok=False)
        with open(os.path.join(S, "gate_G_P1.json"), "w") as fh: json.dump(gate, fh, indent=1)
        raise GateFailure(f"G-P1 failed or not evaluable -- see {S}/gate_G_P1.json")

    # ---- the 240-child cross-release sample
    todo, sample = [], []
    for f in sorted(glob.glob(f"{SAMPLE}/sample_R*.csv")):
        s = pd.read_csv(f)
        for _, r in s.iterrows():
            sample.append((r.release, r.pid))
            if not os.path.exists(os.path.join(d, f"{r.release}_{r.pid}.json")) and \
               not os.path.exists(os.path.join(d, f"{r.release}_{r.pid}.failed.json")):
                todo.append((r.release, r.pid))
    print(f"{len(sample)} sampled children; {len(sample) - len(todo)} done; {len(todo)} to process")
    if MAX_NEW is not None: todo = todo[:MAX_NEW]
    for k, (rel, pid) in enumerate(todo, 1):
        t0 = time.time()
        try:
            fz = CC[(CC.release == rel) & (CC.participant_id == pid)]
            res = run_child(rel, pid, fz); res.pop("bads")
            with open(os.path.join(d, f"{rel}_{pid}.json.tmp"), "w") as fh: json.dump(res, fh)
            os.replace(os.path.join(d, f"{rel}_{pid}.json.tmp"), os.path.join(d, f"{rel}_{pid}.json"))
            msg = f"OK frozen {len(fz)}, caught {sum(r['caught'] for r in res['rows'])}" + (" (RANSAC)" if res["used_ransac"] else "")
        except Exception as e:
            with open(os.path.join(d, f"{rel}_{pid}.failed.json"), "w") as fh:
                json.dump(dict(release=rel, participant_id=pid, reason=f"{type(e).__name__}: {e}"), fh)
            msg = f"FAIL {type(e).__name__}: {e}"
        print(f"[{k}/{len(todo)}] {rel} {pid} {msg} ({time.time() - t0:.0f} s)")

    # ---- summary (new folder)
    S = os.path.join(RUN, f"summary_{time.strftime('%Y%m%d_%H%M%S')}"); os.makedirs(S, exist_ok=False)
    with open(os.path.join(S, "gate_G_P1.json"), "w") as fh: json.dump(gate, fh, indent=1)
    rows = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".json") and not f.endswith(".failed.json"):
            with open(os.path.join(d, f)) as fh: j = json.load(fh)
            for r in j["rows"]:
                rows.append(dict(release=j["release"], participant_id=j["participant_id"], used_ransac=j["used_ransac"], **r))
    D = pd.DataFrame(rows)
    D = D.merge(CC[["release", "participant_id", "channel", "band", "values", "whole_recording", "location"]],
                on=["release", "participant_id", "channel"], how="left")
    D["type"] = D["values"].map(vtype)
    D.to_csv(os.path.join(S, "detector_xrel_per_channel.csv"), index=False)
    fails = [f for f in os.listdir(d) if f.endswith(".failed.json")]

    def block(dim, frame, label):
        out = []
        for key, g in frame.groupby(dim):
            k, n = int(g.caught.sum()), len(g); lo, hi = wilson(k, n)
            out.append(dict(source=label, dimension=dim, level=key, n_channels=n, caught=k,
                            pct=100 * k / n, ci=f"[{100 * lo:.1f}, {100 * hi:.1f}]",
                            caught_without_ransac=int(g.caught_without_ransac.sum())))
        return out
    T = []
    if len(D):
        k, n = int(D.caught.sum()), len(D); lo, hi = wilson(k, n)
        T.append(dict(source="cross-release sample", dimension="all", level="all", n_channels=n, caught=k,
                      pct=100 * k / n, ci=f"[{100 * lo:.1f}, {100 * hi:.1f}]",
                      caught_without_ransac=int(D.caught_without_ransac.sum())))
        for dim in ("band", "type", "release", "location"):
            T += block(dim, D, "cross-release sample")
    er = os.path.join(EARLIER, "detector_stuck_channels.csv")
    if os.path.exists(er):                                       # Release 3, earlier run, labelled as such
        E = pd.read_csv(er)
        E = E.merge(CC[CC.release == "R3"][["participant_id", "channel", "band"]], on=["participant_id", "channel"], how="left")
        T += block("band", E, "Release 3 (earlier run)")
    T = pd.DataFrame(T); T.to_csv(os.path.join(S, "detector_xrel_summary.csv"), index=False)
    crit = D.criteria.str.split(";").explode().replace("", np.nan).dropna().value_counts() if len(D) else pd.Series(dtype=int)
    missed = D[~D.caught] if len(D) else D

    pd.set_option("display.width", 220)
    print(f"\nSUMMARY FOLDER: {S}")
    print(f"GATE G-P1 (gate_G_P1.json): identical = {gate['all_identical']} ({gate['children']} children)")
    print(f"failed children: {len(fails)}")
    print("\nDETECTION (detector_xrel_summary.csv):")
    print(T.round(1).to_string(index=False))
    print("\ncriteria that caught frozen channels (detector_xrel_per_channel.csv): " +
          ", ".join(f"{k} ({v})" for k, v in crit.items()))
    print(f"missed channels: {len(missed)}" + ("" if not len(missed) else "\n" +
          missed[["release", "participant_id", "channel", "band", "type", "values"]].to_string(index=False)))


if __name__ == "__main__":
    main()
