# =============================================================================
# HBN-EEG CROSS-RELEASE STUCK-CHANNEL SCAN        (Colab cell, QC only)
#
# Question: is the undocumented stuck value 193893.78125 (found in 104 of 181
# Release 3 resting recordings; E67 in 89) specific to Release 3, or HBN-wide?
#
# Design (fixed before scanning):
#   * Random sample of N_PER_RELEASE children per release (seed SEED), saved to Drive
#     BEFORE any child is scanned, so the sample can never change on a rerun.
#   * QC only: per channel, non-finite samples, longest run of identical consecutive
#     samples, and that run's value. No filtering, no spectral fitting.
#   * Stuck = identical run >= 1 s (same rule as pipeline v3.1); Cz (reference) exempt.
#   * VALIDATION: first re-scans VALIDATE_N Release 3 children already processed by
#     pipeline v3.1 and requires identical stuck verdicts.
#   * Resumable + crash-safe: one CSV per child, written atomically.
# Run after: drive.mount('/content/drive')
# =============================================================================
import os, json, time, shutil, tempfile, gc
import numpy as np, pandas as pd, requests, scipy.io

S3_ROOT  = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG"
RELEASES = ["R1", "R2", "R4", "R5", "R6", "R7", "R8", "R9"]   # R3 already fully processed
N_PER_RELEASE = 30          # children sampled per release (None = all eligible)
SEED     = 20261021         # fixed sampling seed
OUT_SCAN = "/content/drive/MyDrive/hbn/stuck_scan"
R3_V31   = "/content/drive/MyDrive/hbn/R3_v31"    # for the validation step
VALIDATE_N = 5
INCLUDE_STATUS = ("available", "caution")
KNOWN_CODE = 193893.78125
REF_NAME = "Cz"


def atomic_csv(df, path):
    df.to_csv(path + ".tmp", index=False); os.replace(path + ".tmp", path)


def fetch(url, dest, retries=3):
    last = None
    for k in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                if r.status_code == 404:
                    raise FileNotFoundError(f"404: {url}")
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                enc = r.headers.get("content-encoding", "").lower() not in ("", "identity")
                done = 0
                with open(dest + ".part", "wb") as f:
                    for b in r.iter_content(1 << 20):
                        f.write(b); done += len(b)
            if total and not enc and done != total:
                raise IOError(f"incomplete {done}/{total}")
            os.replace(dest + ".part", dest); return
        except FileNotFoundError:
            raise
        except Exception as e:
            last = e; time.sleep(2 * k)
    raise IOError(f"failed after {retries} attempts: {last}")


def longest_identical_run(x):
    """Identical to pipeline v3.1: longest run of identical consecutive samples."""
    same = np.diff(x) == 0
    if not same.any():
        return 1, float(x[0])
    d = np.diff(np.concatenate([[0], same.astype(np.int8), [0]]))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    k = int(np.argmax(en - st))
    return int(en[k] - st[k] + 1), float(x[st[k]])


def load_any_set(set_path, base_url):
    """Handles top-level-field .set (as in R3), nested EEG struct, and separate .fdt."""
    m = scipy.io.loadmat(set_path, squeeze_me=True, struct_as_record=False)
    if "data" in m:
        data, sf, nb, pn = m["data"], float(m["srate"]), m.get("nbchan"), m.get("pnts")
    elif "EEG" in m:
        E = m["EEG"]; data, sf, nb, pn = E.data, float(E.srate), E.nbchan, E.pnts
    else:
        raise RuntimeError("unsupported .set layout")
    if isinstance(data, str):                                    # separate .fdt
        fdt = os.path.join(os.path.dirname(set_path), data)
        fetch(base_url.rsplit("/", 1)[0] + "/" + data, fdt)
        data = np.fromfile(fdt, dtype=np.float32).reshape(int(pn), int(nb)).T
    return np.asarray(data, dtype=np.float32), sf


def scan_child(release, pid):
    base = f"{S3_ROOT}/cmi_bids_{release}/sub-{pid}/eeg/sub-{pid}_task-RestingState"
    tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
    try:
        sp, cp = os.path.join(tmp, "eeg.set"), os.path.join(tmp, "channels.tsv")
        fetch(base + "_eeg.set", sp); fetch(base + "_channels.tsv", cp)
        X, sf = load_any_set(sp, base + "_eeg.set")
        names = pd.read_csv(cp, sep="\t")["name"].tolist()
        if len(names) != X.shape[0]:
            raise RuntimeError(f"{len(names)} names vs {X.shape[0]} rows")
        n, rows = X.shape[1], []
        for i, c in enumerate(names):
            x = X[i]; L, v = longest_identical_run(x)
            rows.append(dict(release=release, participant_id=pid, channel=c, is_reference=(c == REF_NAME),
                             sfreq=sf, n_samples=n, n_nonfinite=int((~np.isfinite(x)).sum()),
                             longest_run_s=L / sf, run_value=v,
                             run_value_is_known_code=bool(np.isclose(v, KNOWN_CODE, atol=0.05)),
                             whole_recording=bool(L == n),
                             frac_samples_at_known_code=float(np.mean(np.abs(x.astype(np.float64) - KNOWN_CODE) <= 0.05)),
                             stuck=bool(L / sf >= 1.0 and c != REF_NAME)))
        return pd.DataFrame(rows)
    finally:
        shutil.rmtree(tmp, ignore_errors=True); gc.collect()


def get_sample(release, rng_seed):
    """Fixed random sample, saved once; later runs reuse the saved list."""
    path = os.path.join(OUT_SCAN, f"sample_{release}.csv")
    if os.path.exists(path):
        return pd.read_csv(path)["pid"].tolist()
    tsv = os.path.join(OUT_SCAN, f"participants_{release}.tsv")
    fetch(f"{S3_ROOT}/cmi_bids_{release}/participants.tsv", tsv)
    p = pd.read_csv(tsv, sep="\t")
    p = p[p["RestingState"].isin(INCLUDE_STATUS)]
    pids = sorted(p["participant_id"].str.replace("^sub-", "", regex=True))
    rng = np.random.default_rng(rng_seed)
    take = pids if N_PER_RELEASE is None or N_PER_RELEASE >= len(pids) else \
        sorted(rng.choice(pids, N_PER_RELEASE, replace=False).tolist())
    atomic_csv(pd.DataFrame({"pid": take, "release": release, "seed": rng_seed,
                             "n_eligible": len(pids)}), path)
    return take


def run_release(release, pids):
    d = os.path.join(OUT_SCAN, "per_participant"); os.makedirs(d, exist_ok=True)
    todo = [p for p in pids if not os.path.exists(os.path.join(d, f"{release}_{p}.csv"))
            and not os.path.exists(os.path.join(d, f"{release}_{p}.failed.json"))]
    print(f"  {len(pids) - len(todo)} already done, {len(todo)} to scan")
    for k, pid in enumerate(todo, 1):
        t0 = time.time()
        try:
            df = scan_child(release, pid)
            atomic_csv(df, os.path.join(d, f"{release}_{pid}.csv"))
            s = df[df.stuck]
            msg = f"OK  stuck={len(s)} (known code {int(s.run_value_is_known_code.sum())})"
        except Exception as e:
            with open(os.path.join(d, f"{release}_{pid}.failed.json"), "w") as fh:
                json.dump({"release": release, "pid": pid, "reason": f"{type(e).__name__}: {e}"}, fh)
            msg = f"FAIL {type(e).__name__}: {e}"
        print(f"  [{release} {k}/{len(todo)}] {pid}  {msg}  ({time.time() - t0:.0f} s)")


def validate_against_r3():
    """Scanner verdicts must equal pipeline v3.1's stuck verdicts for the same children."""
    C = pd.read_csv(os.path.join(R3_V31, "R3_v31_channels.csv"),
                    usecols=["participant_id", "channel", "invalid_reasons"])
    pids = sorted(C.participant_id.unique())[:VALIDATE_N]
    run_release("R3", pids)
    S = pd.concat([pd.read_csv(os.path.join(OUT_SCAN, "per_participant", f"R3_{p}.csv")) for p in pids])
    C = C[C.participant_id.isin(pids)].copy()
    C["v31_stuck"] = C.invalid_reasons.fillna("").str.contains("stuck")
    M = S.merge(C, on=["participant_id", "channel"])
    agree = (M.stuck == M.v31_stuck).mean()
    print(f"\nVALIDATION vs pipeline v3.1 (R3, {len(pids)} children, {len(M)} channels): "
          f"agreement {agree * 100:.2f}%")
    if agree < 1.0:
        raise RuntimeError("scanner disagrees with pipeline v3.1 -- stop and investigate")


def summarise():
    d = os.path.join(OUT_SCAN, "per_participant")
    files = [f for f in os.listdir(d) if f.endswith(".csv") and not f.startswith("R3_")]
    if not files:
        print("nothing scanned yet"); return
    A = pd.concat([pd.read_csv(os.path.join(d, f)) for f in files], ignore_index=True)
    rows = []
    for rel, g in A.groupby("release"):
        per = g.groupby("participant_id").agg(any_stuck=("stuck", "any"),
              any_code=("run_value_is_known_code", lambda s: bool((s & g.loc[s.index, "stuck"]).any())))
        st = g[g.stuck]
        rows.append(dict(release=rel, children=len(per),
                         with_any_stuck=int(per.any_stuck.sum()),
                         with_known_code=int(per.any_code.sum()),
                         pct_known_code=round(100 * per.any_code.mean(), 1),
                         E67_known_code=int(((g.channel == "E67") & g.stuck & g.run_value_is_known_code).sum()),
                         whole_recording_stuck=int((st.whole_recording).sum()),
                         episodic_stuck=int((~st.whole_recording).sum()),
                         other_values=";".join(f"{v:.5g}" for v in sorted(st[~st.run_value_is_known_code].run_value.unique())[:5])))
    T = pd.DataFrame(rows)
    atomic_csv(T, os.path.join(OUT_SCAN, "summary_by_release.csv"))
    st = A[A.stuck & A.run_value_is_known_code]
    top = st.channel.value_counts().head(8)
    fails = [f for f in os.listdir(d) if f.endswith(".failed.json") and not f.startswith("R3_")]
    pd.set_option("display.width", 200)
    print("\nSTUCK-CHANNEL PREVALENCE BY RELEASE"); print(T.to_string(index=False))
    tot = A.groupby(["release", "participant_id"]).apply(
        lambda g: bool((g.stuck & g.run_value_is_known_code).any()), include_groups=False)
    print(f"\nALL SCANNED RELEASES: {int(tot.sum())} of {len(tot)} children "
          f"({100 * tot.mean():.1f}%) have >=1 channel stuck at {KNOWN_CODE}")
    print("most often stuck at that value: " + ", ".join(f"{c} ({k})" for c, k in top.items()))
    print(f"failed downloads/reads: {len(fails)}")


def main():
    if OUT_SCAN.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    os.makedirs(os.path.join(OUT_SCAN, "per_participant"), exist_ok=True)
    validate_against_r3()
    for k, rel in enumerate(RELEASES):
        try:
            pids = get_sample(rel, SEED + k)
        except FileNotFoundError:
            print(f"\n{rel}: participants.tsv not found at cmi_bids_{rel} -- release skipped"); continue
        print(f"\n{rel}: {len(pids)} children in the fixed sample")
        run_release(rel, pids)
    summarise()


if __name__ == "__main__":
    main()
