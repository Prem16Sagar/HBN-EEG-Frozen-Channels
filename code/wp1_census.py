# =============================================================================
# WP1 -- CENSUS OF FROZEN CHANNELS, ALL HBN-EEG RELEASES     (Colab cell)
# Implements revision_B_spec.md, WP1 (approved).
#
# * Every RestingState recording ("available"/"caution") in cmi_bids_R1..R11.
#   A release whose participants.tsv is not found is reported as ABSENT.
# * Per channel, on the full continuous recording (stored values, no filtering):
#   every maximal run of identical consecutive samples lasting >= 1 s (the v3.1 rule).
# * Resumable: each child's result is written ONCE to RUN/per_child/ and never rewritten.
#   Every summary goes to a NEW folder RUN/summary_<timestamp>/. Nothing is overwritten.
# * Gates: G-C1 (Release 3 vs pipeline v3.1) and G-C2 (vs the earlier 240-child scan)
#   must pass before any summary is written.
# Run after: drive.mount('/content/drive')
# =============================================================================
import os, json, time, shutil, tempfile, gc
import numpy as np, pandas as pd, requests, scipy.io
from scipy import stats

DRIVE    = "/content/drive/MyDrive/hbn"
RUN      = f"{DRIVE}/revision_B/census/census_v1"
S3_ROOT  = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG"
RELEASES = [f"R{i}" for i in range(1, 12)]
V31_CHANNELS = f"{DRIVE}/R3_v31/R3_v31_channels.csv"
EARLIER_SCAN = f"{DRIVE}/stuck_scan/per_participant"
INCLUDE_STATUS = ("available", "caution")
MAX_NEW  = None                       # children to scan this session (None = all remaining)
KNOWN    = 193893.78125
REF_NAME = "Cz"
PERIPH = {"E48","E119","E127","E126","E113","E49","E56","E107","E99","E63",
          "E94","E68","E43","E120","E88","E73","E81","E125","E128"}


class GateFailure(RuntimeError):
    pass


def fetch(url, dest, retries=3):
    last = None
    for k in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                if r.status_code == 404:
                    raise FileNotFoundError(f"404: {url}")
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
    """Top-level-field .set (as in R3), nested EEG struct, or separate .fdt."""
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
    return np.asarray(data, dtype=np.float32), sf


def runs_ge_1s(x, sf):
    """All maximal runs of identical consecutive samples with duration >= 1 s.
    Run length counted in samples exactly as pipeline v3.1 (en - st + 1)."""
    same = np.diff(x) == 0
    if not same.any():
        return 1, []
    d = np.diff(np.concatenate([[0], same.astype(np.int8), [0]]))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    L = en - st + 1
    longest = int(L.max())
    keep = np.where(L / sf >= 1.0)[0]
    return longest, [dict(start_s=float(st[k] / sf), duration_s=float(L[k] / sf), n_samples=int(L[k]),
                          value=float(x[st[k]])) for k in keep]


def band(longest_s, whole):
    if whole: return "whole recording"
    if longest_s >= 30: return "[30 s, whole)"
    if longest_s >= 5: return "[5, 30) s"
    return "[1, 5) s"


def scan_child(rel, pid):
    base = f"{S3_ROOT}/cmi_bids_{rel}/sub-{pid}/eeg/sub-{pid}_task-RestingState"
    tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
    try:
        sp, cp = os.path.join(tmp, "eeg.set"), os.path.join(tmp, "channels.tsv")
        fetch(base + "_eeg.set", sp); fetch(base + "_channels.tsv", cp)
        X, sf = load_any_set(sp, base + "_eeg.set")
        names = pd.read_csv(cp, sep="\t")["name"].tolist()
        if len(names) != X.shape[0]:
            raise RuntimeError(f"{len(names)} channel names vs {X.shape[0]} data rows")
        n = X.shape[1]; chans, longest_all = [], {}
        for i, c in enumerate(names):
            x = X[i]; L, eps = runs_ge_1s(x, sf); longest_all[c] = L / sf
            if c == REF_NAME or not eps:
                continue
            whole = len(eps) == 1 and eps[0]["n_samples"] == n
            ls = max(e["duration_s"] for e in eps)
            chans.append(dict(channel=c, location="peripheral" if c in PERIPH else "scalp",
                              n_episodes=len(eps), total_frozen_s=float(sum(e["duration_s"] for e in eps)),
                              longest_s=ls, whole_recording=bool(whole), band=band(ls, whole),
                              values=sorted({repr(e["value"]) for e in eps}),
                              n_nonfinite=int((~np.isfinite(x)).sum()), episodes=eps))
        return dict(release=rel, participant_id=pid, sfreq=sf, n_samples=int(n), duration_s=n / sf,
                    n_channels=len(names), frozen=chans, longest_run_s_all=longest_all)
    finally:
        shutil.rmtree(tmp, ignore_errors=True); gc.collect()


def atomic_json(obj, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh)
    os.replace(path + ".tmp", path)


def participants(rel):
    d = os.path.join(RUN, "participants"); os.makedirs(d, exist_ok=True)
    p = os.path.join(d, f"participants_{rel}.tsv")
    if not os.path.exists(p):
        fetch(f"{S3_ROOT}/cmi_bids_{rel}/participants.tsv", p)
    t = pd.read_csv(p, sep="\t")
    t = t[t["RestingState"].isin(INCLUDE_STATUS)]
    return sorted(t["participant_id"].str.replace("^sub-", "", regex=True))


def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


# ----------------------------------------------------------------------------- gates
def gates(J):
    out = {}
    # G-C1: Release 3 vs pipeline v3.1, every channel
    if os.path.exists(V31_CHANNELS):
        V = pd.read_csv(V31_CHANNELS, usecols=["participant_id", "channel", "invalid_reasons", "longest_identical_run_s"])
        V["v31_frozen"] = V.invalid_reasons.fillna("").str.contains("stuck")
        rows = []
        for j in J:
            if j["release"] != "R3": continue
            fz = {c["channel"] for c in j["frozen"]}
            for c, ls in j["longest_run_s_all"].items():
                rows.append(dict(participant_id=j["participant_id"], channel=c, census_frozen=c in fz,
                                 census_longest=ls, sf=j["sfreq"]))
        M = pd.DataFrame(rows).merge(V, on=["participant_id", "channel"], how="inner") if rows else pd.DataFrame()
        if len(M):
            vmis = int((M.census_frozen != M.v31_frozen).sum())
            dmis = int((np.abs(M.census_longest - M.longest_identical_run_s) > 1.0 / M.sf + 1e-9).sum())
            out["G-C1"] = dict(children=int(M.participant_id.nunique()), channels=len(M),
                               verdict_mismatches=vmis, duration_mismatches=dmis, passed=vmis == 0 and dmis == 0)
        else:
            out["G-C1"] = dict(children=0, passed=None, note="no Release 3 children scanned yet")
    # G-C2: the 240 earlier-sample children vs the earlier scan files
    if os.path.isdir(EARLIER_SCAN):
        rows = []
        for j in J:
            f = os.path.join(EARLIER_SCAN, f"{j['release']}_{j['participant_id']}.csv")
            if j["release"] == "R3" or not os.path.exists(f): continue
            E = pd.read_csv(f, usecols=["channel", "stuck", "longest_run_s"])
            fz = {c["channel"] for c in j["frozen"]}
            for _, r in E.iterrows():
                cl = j["longest_run_s_all"].get(r.channel, np.nan)
                rows.append(dict(cv=r.channel in fz, ev=str(r.stuck).lower() == "true",
                                 dd=abs(cl - r.longest_run_s) > 1.0 / j["sfreq"] + 1e-9, pid=j["participant_id"]))
        if rows:
            M = pd.DataFrame(rows)
            out["G-C2"] = dict(children=int(M.pid.nunique()), channels=len(M),
                               verdict_mismatches=int((M.cv != M.ev).sum()), duration_mismatches=int(M.dd.sum()),
                               passed=bool((M.cv == M.ev).all() and not M.dd.any()))
        else:
            out["G-C2"] = dict(children=0, passed=None, note="no earlier-sample children scanned yet")
    return out


# ----------------------------------------------------------------------------- summary
def summarise():
    d = os.path.join(RUN, "per_child")
    J = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".json") and not f.endswith(".failed.json"):
            with open(os.path.join(d, f)) as fh: J.append(json.load(fh))
    fails = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".failed.json"):
            with open(os.path.join(d, f)) as fh: fails.append(json.load(fh))
    G = gates(J)
    S = os.path.join(RUN, f"summary_{time.strftime('%Y%m%d_%H%M%S')}"); os.makedirs(S, exist_ok=False)
    atomic_json(G, os.path.join(S, "gates.json"))
    if any(v.get("passed") is False for v in G.values()):
        raise GateFailure(f"census gate failed -- see {S}/gates.json: {json.dumps(G)}")

    ch, ep, rec = [], [], []
    for j in J:
        fz = j["frozen"]
        for c in fz:
            ch.append(dict(release=j["release"], participant_id=j["participant_id"], channel=c["channel"],
                           location=c["location"], n_episodes=c["n_episodes"], total_frozen_s=c["total_frozen_s"],
                           longest_s=c["longest_s"], whole_recording=c["whole_recording"], band=c["band"],
                           values=";".join(c["values"]), n_nonfinite=c["n_nonfinite"],
                           at_known_code=any(float(v) == KNOWN for v in c["values"])))
            for e in c["episodes"]:
                ep.append(dict(release=j["release"], participant_id=j["participant_id"], channel=c["channel"],
                               start_s=e["start_s"], duration_s=e["duration_s"], value=repr(e["value"])))
        rec.append(dict(release=j["release"], participant_id=j["participant_id"], duration_s=j["duration_s"],
                        n_frozen=len(fz), n_frozen_scalp=sum(c["location"] == "scalp" for c in fz),
                        n_whole=sum(c["whole_recording"] for c in fz), n_episodic=sum(not c["whole_recording"] for c in fz),
                        any_known_code=any(float(v) == KNOWN for c in fz for v in c["values"]),
                        any_frozen=len(fz) > 0, any_episodic=any(not c["whole_recording"] for c in fz)))
    CH, EP, REC = pd.DataFrame(ch), pd.DataFrame(ep), pd.DataFrame(rec)
    CH.to_csv(os.path.join(S, "census_channels.csv"), index=False)
    EP.to_csv(os.path.join(S, "census_episodes.csv"), index=False)
    REC.to_csv(os.path.join(S, "census_recordings.csv"), index=False)
    pd.DataFrame(fails).to_csv(os.path.join(S, "census_failures.csv"), index=False)

    status = json.load(open(os.path.join(RUN, "release_status.json")))
    rows = []
    for rel in RELEASES:
        g = REC[REC.release == rel] if len(REC) else REC; n = len(g)
        row = dict(release=rel, status=status.get(rel, {}).get("status", "not attempted"),
                   listed=status.get(rel, {}).get("listed", 0), scanned=n,
                   failed=sum(f["release"] == rel for f in fails))
        for key, col in (("known_code", "any_known_code"), ("any_frozen", "any_frozen"), ("episodic", "any_episodic")):
            k = int(g[col].sum()) if n else 0; lo, hi = wilson(k, n)
            row.update({f"n_{key}": k, f"pct_{key}": 100 * k / n if n else np.nan,
                        f"ci_{key}": f"{100 * lo:.1f}-{100 * hi:.1f}" if n else ""})
        cc = CH[(CH.release == rel) & CH.at_known_code] if len(CH) else CH
        row["dominant_channels"] = ", ".join(f"{c} ({k})" for c, k in
                                             cc.groupby("channel").participant_id.nunique().sort_values(ascending=False).head(3).items()) if len(cc) else ""
        rows.append(row)
    N = len(REC); tot = dict(release="ALL (census)", status="", listed=sum(r["listed"] for r in rows), scanned=N,
                             failed=len(fails), dominant_channels="")
    for key, col in (("known_code", "any_known_code"), ("any_frozen", "any_frozen"), ("episodic", "any_episodic")):
        k = int(REC[col].sum()) if N else 0; lo, hi = wilson(k, N)
        tot.update({f"n_{key}": k, f"pct_{key}": 100 * k / N if N else np.nan, f"ci_{key}": f"{100 * lo:.1f}-{100 * hi:.1f}" if N else ""})
    T = pd.DataFrame(rows + [tot]); T.to_csv(os.path.join(S, "census_summary.csv"), index=False)

    extra = {}
    scanned = [r for r in rows if r["scanned"] > 0]
    for key in ("known_code", "any_frozen"):
        tab = np.array([[r[f"n_{key}"] for r in scanned], [r["scanned"] - r[f"n_{key}"] for r in scanned]])
        if tab.shape[1] > 1 and (tab.sum(1) > 0).all():
            chi2, p, dof, _ = stats.chi2_contingency(tab)
            extra[f"chi2_{key}"] = dict(chi2=float(chi2), dof=int(dof), p=float(p))
    if len(EP):
        extra["episode_duration_s"] = dict(n=len(EP), median=float(EP.duration_s.median()),
                                           q25=float(EP.duration_s.quantile(.25)), q75=float(EP.duration_s.quantile(.75)))
        extra["channels_by_band"] = CH.band.value_counts().to_dict()
    atomic_json(extra, os.path.join(S, "census_extra.json"))
    return S, T, G, extra


def main():
    if RUN.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    d = os.path.join(RUN, "per_child"); os.makedirs(d, exist_ok=True)
    sp = os.path.join(RUN, "release_status.json")
    status = json.load(open(sp)) if os.path.exists(sp) else {}
    todo = []
    for rel in RELEASES:
        try:
            pids = participants(rel); status[rel] = dict(status="found", listed=len(pids))
            todo += [(rel, p) for p in pids if not os.path.exists(os.path.join(d, f"{rel}_{p}.json"))
                     and not os.path.exists(os.path.join(d, f"{rel}_{p}.failed.json"))]
        except FileNotFoundError:
            status[rel] = dict(status="ABSENT: participants.tsv not found", listed=0)
    atomic_json(status, sp)
    for rel, s in status.items():
        print(f"  {rel}: {s['status']}, {s['listed']} listed")
    total = sum(s["listed"] for s in status.values())
    print(f"\n{total} recordings listed; {len(todo)} still to scan" + (f" (this session: {MAX_NEW})" if MAX_NEW else ""))
    if MAX_NEW is not None: todo = todo[:MAX_NEW]
    t_start = time.time()
    for k, (rel, pid) in enumerate(todo, 1):
        t0 = time.time()
        try:
            res = scan_child(rel, pid); atomic_json(res, os.path.join(d, f"{rel}_{pid}.json"))
            msg = f"OK frozen={len(res['frozen'])}"
        except Exception as e:
            atomic_json(dict(release=rel, participant_id=pid, reason=f"{type(e).__name__}: {e}"),
                        os.path.join(d, f"{rel}_{pid}.failed.json"))
            msg = f"FAIL {type(e).__name__}: {e}"
        if k % 25 == 0 or k == len(todo) or msg.startswith("FAIL"):
            eta = (time.time() - t_start) / k * (len(todo) - k) / 60
            print(f"  [{k}/{len(todo)}] {rel} {pid} {msg} ({time.time() - t0:.0f} s; ETA {eta:.0f} min)")
    S, T, G, extra = summarise()
    pd.set_option("display.width", 250)
    print(f"\nSUMMARY FOLDER: {S}")
    print("GATES (gates.json):"); [print(f"  {k}: {v}") for k, v in G.items()]
    print("\nCENSUS BY RELEASE (census_summary.csv):")
    print(T[["release", "status", "listed", "scanned", "failed", "n_known_code", "pct_known_code", "ci_known_code",
             "n_any_frozen", "pct_any_frozen", "ci_any_frozen", "n_episodic", "pct_episodic", "ci_episodic",
             "dominant_channels"]].round(1).to_string(index=False))
    print("\nEXTRA (census_extra.json):"); print(json.dumps(extra, indent=1))


if __name__ == "__main__":
    main()
