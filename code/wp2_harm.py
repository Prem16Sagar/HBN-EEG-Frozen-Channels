# =============================================================================
# WP2 -- COHORT-LEVEL HARM OF FROZEN CHANNELS, RELEASE 3       (Colab cell)
# Implements revision_B_spec.md, WP2 (approved).
#
# Processing is pipeline v3.1's own code, imported unchanged from pilot_v31.py:
# load_set, ec_blocks, preprocess_channel (1 Hz high-pass + notches, per channel,
# BEFORE referencing), Welch per eyes-closed block, fit_aperiodic (fooof 1.1.1, 2-40 Hz),
# spectral_features (1-110 Hz RMS, line bins excluded). Quality-control flags are read
# unchanged from R3_v31_channels.csv. ONLY the average reference is varied:
#   correct      valid, non-high-variance scalp + Cz            (v3.1)
#   wrong        ALL scalp channels + Cz                        (reference before detection)
#   frozen_only  correct set + frozen scalp channels            (frozen contribution alone)
# Gate G-H1: the "correct" exponent must reproduce the saved scalp__ap2_exponent within
# 1e-6 for every child; the first mismatch stops the run.
# Resumable: per-child results written once to RUN/per_child/; summaries to new folders.
# Run after: drive.mount('/content/drive'). pilot_v31.py must be in CODE_DIR.
# =============================================================================
import os, sys, json, time, shutil, tempfile, gc, warnings
import numpy as np, pandas as pd, requests
from scipy.signal import welch
import statsmodels.api as sm, patsy

DRIVE    = "/content/drive/MyDrive/hbn"
CODE_DIR = DRIVE                                   # folder holding pilot_v31.py
R3_DIR   = f"{DRIVE}/R3_v31"
RUN      = f"{DRIVE}/revision_B/harm/harm_v1"
S3_R3    = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG/cmi_bids_R3"
MAX_NEW  = None
N_BOOT, SEED, RMIN, TOL_H1 = 2000, 20261023, 0.85, 1e-6
ERR_THRESHOLD = 0.05


class GateFailure(RuntimeError):
    pass


def _bool(s):
    return s if s.dtype == bool else s.astype(str).str.strip().str.lower().isin(("true", "1"))


def fetch(url, dest, retries=3):
    last = None
    for k in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(dest + ".part", "wb") as f:
                    for b in r.iter_content(1 << 20): f.write(b)
            os.replace(dest + ".part", dest); return
        except Exception as e:
            last = e; time.sleep(2 * k)
    raise IOError(f"failed after {retries} attempts: {last}")


def atomic_json(obj, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, default=float)
    os.replace(path + ".tmp", path)


def process_child(v, pid, F, saved_S, folder):
    """F: this child's rows of R3_v31_channels.csv (flags used unchanged)."""
    from scipy.signal import butter
    X, sf, names, evp = v.load_set(folder, pid)
    n = X.shape[1]
    blocks = v.ec_blocks(evp, sf, n)
    sos_hp = butter(4, 1.0, btype="highpass", fs=sf, output="sos")
    F = F.set_index("channel").loc[names]                          # align to file channel order
    Y = np.zeros((len(names), n))
    for i in range(len(names)):                                    # exactly as v3.1: finite channels only
        x = X[i]
        if np.isfinite(x).all():
            Y[i] = v.preprocess_channel(x.astype(np.float64), sf, sos_hp)
    loc = F.location.values; inval = F.is_invalid.values; inref = F.in_reference.values
    frozen = (loc == "scalp") & F.invalid_reasons.fillna("").str.contains("stuck").values
    cz = names.index(v.REF_NAME)
    sets = {"correct": np.where(inref)[0],
            "wrong": np.array(sorted(set(np.where(loc == "scalp")[0]) | {cz})),
            "frozen_only": np.array(sorted(set(np.where(inref)[0]) | set(np.where(frozen)[0])))}
    summ = np.where((loc == "scalp") & ~inval & inref)[0]          # v3.1 summary set (unchanged)
    valid_scalp = np.where((loc == "scalp") & ~inval)[0]           # channels used for the RMS ratio
    out = {}
    for cond, idx in sets.items():
        Yc = Y - Y[idx].mean(axis=0, keepdims=True)
        Pb = []
        for s, e in blocks:
            f, P = welch(Yc[:, s:e], fs=sf, nperseg=int(2 * sf), noverlap=int(sf), window="hann", axis=-1)
            Pb.append(P)
        Pm = np.mean(Pb, axis=0)
        fits = [v.fit_aperiodic(f, Pm[i], 2) for i in summ]
        ex = pd.Series([d["exponent"] for d in fits]); r2 = pd.Series([d["r2"] for d in fits])
        ok = r2 >= RMIN
        rms = pd.Series([v.spectral_features(f, Pm[i])["rms_1_110_uV"] for i in valid_scalp])
        out[cond] = dict(exponent=float(ex[ok].median()) if ok.any() else np.nan,
                         n_fit_fail=int((~ok).sum()), median_rms_1_110_uV=float(rms.median()),
                         n_reference=int(len(idx)))
    if not abs(out["correct"]["exponent"] - saved_S) <= TOL_H1:
        raise GateFailure(f"G-H1: {pid} correct-reference exponent {out['correct']['exponent']!r} "
                          f"differs from saved {saved_S!r}")
    longest = F.longest_identical_run_s.values
    whole = frozen & (longest >= n / sf - 1.0 / sf)
    n_whole, n_epi = int(whole.sum()), int((frozen & ~whole).sum())
    group = "no frozen scalp channel" if frozen.sum() == 0 else (
        ">=1 episodic frozen scalp channel" if n_epi else "whole-recording frozen scalp channels only")
    frozen_periph = int(((loc == "peripheral") & F.invalid_reasons.fillna("").str.contains("stuck").values).sum())
    res = dict(participant_id=pid, group=group, n_frozen_scalp_whole=n_whole, n_frozen_scalp_episodic=n_epi,
               n_frozen_peripheral=frozen_periph, saved_exponent=saved_S, gate_abs_diff=abs(out["correct"]["exponent"] - saved_S))
    for cond in ("correct", "wrong", "frozen_only"):
        for k, val in out[cond].items():
            res[f"{cond}__{k}"] = val
    for cond in ("wrong", "frozen_only"):
        res[f"{cond}__inflation_ratio"] = out[cond]["median_rms_1_110_uV"] / out["correct"]["median_rms_1_110_uV"]
        res[f"{cond}__exponent_error"] = out[cond]["exponent"] - out["correct"]["exponent"]
    return res


def summarise():
    d = os.path.join(RUN, "per_child")
    def _read(f):
        with open(os.path.join(d, f)) as fh:
            return json.load(fh)
    J = [_read(f) for f in sorted(os.listdir(d)) if f.endswith(".json") and not f.endswith(".failed.json")]
    fails = [_read(f) for f in sorted(os.listdir(d)) if f.endswith(".failed.json")]
    S = os.path.join(RUN, f"summary_{time.strftime('%Y%m%d_%H%M%S')}"); os.makedirs(S, exist_ok=False)
    H = pd.DataFrame(J)
    A = pd.read_csv(f"{R3_DIR}/analysis/analysis_sample.csv", usecols=["participant_id", "age", "sex"])
    H = H.merge(A, on="participant_id", how="left")
    H.to_csv(os.path.join(S, "harm_per_child.csv"), index=False)
    pd.DataFrame(fails).to_csv(os.path.join(S, "harm_failures.csv"), index=False)
    rng = np.random.default_rng(SEED)

    def med_ci(x):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        if len(x) < 2: return (np.nan, np.nan)
        b = [np.median(x[rng.integers(0, len(x), len(x))]) for _ in range(N_BOOT)]
        return tuple(np.percentile(b, [2.5, 97.5]))

    rows = []
    groups = ["ALL"] + sorted(H.group.unique())
    for grp in groups:
        g = H if grp == "ALL" else H[H.group == grp]
        for cond in ("wrong", "frozen_only"):
            ir, er = g[f"{cond}__inflation_ratio"], g[f"{cond}__exponent_error"]
            lo_i, hi_i = med_ci(ir); lo_e, hi_e = med_ci(er)
            rows.append(dict(group=grp, reference=cond, n=len(g),
                             inflation_median=float(ir.median()), inflation_q25=float(ir.quantile(.25)),
                             inflation_q75=float(ir.quantile(.75)), inflation_ci=f"[{lo_i:.4g}, {hi_i:.4g}]",
                             inflation_max=float(ir.max()),
                             error_median=float(er.median()), error_q25=float(er.quantile(.25)),
                             error_q75=float(er.quantile(.75)), error_ci=f"[{lo_e:.4g}, {hi_e:.4g}]",
                             pct_abs_error_gt_0_05=float(100 * (er.abs() > ERR_THRESHOLD).mean()),
                             extra_fit_failures_median=float((g[f"{cond}__n_fit_fail"] - g["correct__n_fit_fail"]).median())))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(S, "harm_summary.csv"), index=False)

    # age slope, wrong vs correct (Huber; coefficient by name; paired bootstrap)
    q = H[["correct__exponent", "wrong__exponent", "age", "sex"]].dropna()
    X = patsy.dmatrix("1 + age + C(sex)", q, return_type="dataframe", eval_env=patsy.EvalEnvironment([{}]))
    j = list(X.columns).index("age"); Xv = X.values
    yc, yw = q["correct__exponent"].values, q["wrong__exponent"].values
    def hub(y, Xm):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); return sm.RLM(y, Xm, M=sm.robust.norms.HuberT()).fit().params[j]
    bc, bw = float(hub(yc, Xv)), float(hub(yw, Xv)); B = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(yc), len(yc))
        try: B.append(hub(yw[i], Xv[i]) - hub(yc[i], Xv[i]))
        except Exception: pass
    lo, hi = np.percentile(B, [2.5, 97.5])
    R = dict(gate_G_H1=dict(children=len(H), max_abs_diff=float(H.gate_abs_diff.max()), tolerance=TOL_H1),
             failures=len(fails), age_slope=dict(n=len(q), correct=bc, wrong=bw, difference=bw - bc,
                                                  diff_ci=[float(lo), float(hi)], boot_resamples=len(B)),
             seed=SEED, n_boot=N_BOOT, error_threshold=ERR_THRESHOLD)
    atomic_json(R, os.path.join(S, "harm_results.json"))
    return S, T, R


def main():
    if RUN.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    sys.path.insert(0, CODE_DIR)
    import pilot_v31 as v
    d = os.path.join(RUN, "per_child"); os.makedirs(d, exist_ok=True)
    A = pd.read_csv(f"{R3_DIR}/analysis/analysis_sample.csv", usecols=["participant_id", "scalp__ap2_exponent"])
    CH = pd.read_csv(f"{R3_DIR}/R3_v31_channels.csv",
                     usecols=["participant_id", "channel", "location", "is_invalid", "in_reference",
                              "invalid_reasons", "longest_identical_run_s"])
    CH["is_invalid"], CH["in_reference"] = _bool(CH.is_invalid), _bool(CH.in_reference)
    todo = [p for p in A.participant_id if not os.path.exists(os.path.join(d, f"{p}.json"))
            and not os.path.exists(os.path.join(d, f"{p}.failed.json"))]
    print(f"{len(A)} children; {len(A) - len(todo)} done; {len(todo)} to process")
    if MAX_NEW is not None: todo = todo[:MAX_NEW]
    t_all = time.time()
    for k, pid in enumerate(todo, 1):
        t0 = time.time(); tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
        try:
            base = f"sub-{pid}_task-RestingState"
            for sfx in ("_eeg.set", "_events.tsv", "_channels.tsv"):
                fetch(f"{S3_R3}/sub-{pid}/eeg/{base}{sfx}", os.path.join(tmp, base + sfx))
            saved = float(A.loc[A.participant_id == pid, "scalp__ap2_exponent"].iloc[0])
            res = process_child(v, pid, CH[CH.participant_id == pid], saved, tmp)
            atomic_json(res, os.path.join(d, f"{pid}.json"))
            msg = (f"OK {res['group'][:22]:<22} infl wrong {res['wrong__inflation_ratio']:.2f} "
                   f"err {res['wrong__exponent_error']:+.3f}")
        except GateFailure:
            raise                                                        # method divergence: stop everything
        except Exception as e:
            atomic_json(dict(participant_id=pid, reason=f"{type(e).__name__}: {e}"), os.path.join(d, f"{pid}.failed.json"))
            msg = f"FAIL {type(e).__name__}: {e}"
        finally:
            shutil.rmtree(tmp, ignore_errors=True); gc.collect()
        eta = (time.time() - t_all) / k * (len(todo) - k) / 60
        print(f"[{k}/{len(todo)}] {pid} {msg} ({time.time() - t0:.0f} s, ETA {eta:.0f} min)")
    S, T, R = summarise()
    pd.set_option("display.width", 250)
    print(f"\nSUMMARY FOLDER: {S}")
    print(f"GATE G-H1 (harm_results.json): {R['gate_G_H1']}")
    print("\nHARM BY GROUP (harm_summary.csv):")
    print(T.round(4).to_string(index=False))
    a = R["age_slope"]
    print(f"\nAGE SLOPE (harm_results.json): correct {a['correct']:+.5f}/yr, wrong {a['wrong']:+.5f}/yr, "
          f"difference {a['difference']:+.5f} {np.round(a['diff_ci'], 5).tolist()}, n={a['n']}")


if __name__ == "__main__":
    main()
