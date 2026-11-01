# =============================================================================
# HBN R3 PERIPHERAL-ACTIVITY PIPELINE -- v3.1  (importable module; pasted as a cell it
# runs the 3-participant pilot)
# Feature extraction + QC only. No modelling, no inference, no clinical claims.
# Read-only on source files; writes to OUT only. Run AFTER setup_hbn_files.py.
#
# Changes from v2 (each motivated by a documented diagnostic):
#   * Frozen channel QC (stuck runs, flat/non-finite, notched high-variance, Cz = reference)
#   * Aperiodic exponent + offset (fooof 1.1.1, fixed mode) -- the PRIMARY OUTCOME
#       primary  : 1 Hz high-pass, fit 2-40 Hz  (drift-invariant in simulation)
#       secondary: 1 Hz high-pass, fit 1-40 Hz  (literature-comparable sensitivity arm)
#   * High-frequency noise-floor descriptor (70-110 Hz relative to 150-200 Hz)
#   * Pure numpy/scipy pipeline (no MNE): loads the .set with scipy, as verified
# =============================================================================
import os, sys, json, platform, subprocess, warnings
import numpy as np, pandas as pd, scipy, scipy.io
from scipy.signal import butter, sosfiltfilt, welch, iirnotch, filtfilt

try:                                    # fooof needs only numpy+scipy; --no-deps guarantees
    import fooof                        # no other Colab package is touched
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "--no-deps", "fooof==1.1.1"])
    import fooof
from fooof import FOOOF
warnings.filterwarnings("ignore", category=DeprecationWarning)

PARTICIPANTS = {"NDARAA948VFH": "/content/hbn_one",
                "NDARAC853DTE": "/content/hbn_two",
                "NDARAD774HAZ": "/content/hbn_three"}
OUT = "/content/pilot_out_v31"

GROUPS = {"infraorbital": ["E126", "E127"],
          "face_jaw": ["E48", "E119", "E43", "E120", "E125", "E128", "E49", "E113"],
          "neck": ["E56", "E107", "E63", "E99", "E68", "E94", "E73", "E88", "E81"]}
PERIPH = set(sum(GROUPS.values(), []))
# A-priori frontal scalp group: scalp channels in the anterior 25% of the scalp's
# antero-posterior extent (GSN-HydroCel-129 geometry, y >= 52.6 mm). Defined from geometry
# BEFORE inspecting pilot results; 5 of the 7 pilot HF-elevated channels fall inside it.
FRONTAL = {"E1","E2","E3","E4","E5","E8","E9","E10","E11","E12","E14","E15","E16","E17",
           "E18","E19","E20","E21","E22","E23","E24","E25","E26","E27","E32","E118","E123","E124"}
REF_NAME = "Cz"
KNOWN_STUCK = 193893.78125

PARAMS = {
    "version": "v3.1",
    "v31_changes": ["frontal scalp group (a-priori geometry) and frontal-vs-nonfrontal gradient",
                    "exponent reported with AND without the R^2 filter, plus fit-failure counts "
                    "(the filter censors contaminated channels)",
                    "block-level SD replaces CV (CV is meaningless for slopes near zero)",
                    "knee-mode aperiodic arm, fit 1-40 Hz"],
    "units": "stored values interpreted as microvolts per BIDS channels.tsv; no rescaling; "
             "acquisition-to-export conversion provenance unverified",
    "ec_segmentation": "instructed close->open cue pairs from events.tsv, 2 s trimmed at each edge",
    "highpass": "Butterworth order 4, 1.0 Hz, zero-phase (sosfiltfilt)",
    "notch": "iirnotch Q=30 at 60/120/180/240 Hz, zero-phase (filtfilt), applied after high-pass",
    "qc_stuck": "longest run of identical consecutive samples >= 1 s -> invalid (Cz exempt)",
    "qc_flat": "non-finite samples, or SD < 1e-3 uV on high-passed eyes-closed data -> invalid",
    "qc_high_variance": "robust z of log-SD on notched+high-passed eyes-closed data with the "
                        "ALPHA BAND REMOVED (Butterworth band-stop 7-14 Hz) > +5 (one-sided), among "
                        "valid scalp channels -> excluded from reference and scalp summaries. Alpha "
                        "is removed because strong eyes-closed posterior alpha is physiological, and "
                        "in simulation a broadband rule excluded the whole posterior alpha cluster.",
    "qc_hf_elevated": "robust z of log 65-100 Hz power > +5 among valid scalp -> RECORDED ONLY "
                      "(this activity is the phenomenon under study)",
    "reference": "average of valid, non-high-variance scalp channels plus Cz; peripheral "
                 "channels never in the reference",
    "psd": "Welch, 2 s Hann windows, 50% overlap, per eyes-closed block, mean across blocks",
    "aperiodic_primary": "fooof 1.1.1 fixed mode, fit 2-40 Hz, peak_width_limits [1,8], "
                         "max_n_peaks 6",
    "aperiodic_secondary": "same, fit 1-40 Hz (literature-comparable; sensitivity arm)",
    "aperiodic_choice_rationale": "In simulation with realistic slow electrode drift (random "
        "walk low-passed at 0.1 Hz, 0-30000 uV p2p) and true exponents 1.0-2.5: 1 Hz HP + fit "
        "2-40 gave worst |bias| 0.030 and drift sensitivity 0.001; 1 Hz HP + fit 1-40 gave "
        "0.091 / 0.014; 0.5 Hz HP + fit 1-40 gave 0.078 / 0.042; 0.1 Hz HP or no HP failed "
        "(bias up to 1.7) through drift leakage.",
    "aperiodic_knee_arm": "fooof 1.1.1 knee mode, fit 1-40 Hz (sensitivity: real 2-40 vs 1-40 "
                          "fixed-mode exponents diverged more than simulation predicted)",
    "aperiodic_fit_quality_min_r2": 0.85,
    "focal_hf_metric": "fraction of a group's channels whose rel_70_110 (or hf_frac_30_110) "
                       "exceeds the 90th percentile of non-frontal valid scalp. Sensitive to muscle on "
                       "a minority of channels, which group medians miss. CAUTION: its ~0.10 null holds "
                       "only for exchangeable channel groups (frontal vs non-frontal scalp). Peripheral "
                       "channels sit further from cortex and have flatter spectra even without muscle, "
                       "so for peripheral groups the null must be estimated empirically, not assumed.",
    "notch_exclusion_halfwidth_hz": 2.0,
    "features": {"rel_30_45": "P(30-45)/P(1-45)", "rel_70_110": "P(70-110)/P(1-110)",
                 "hf_frac_30_110": "P(30-110)/P(1-110)",
                 "slope_30_90": "log-log OLS slope 30-90 Hz, 55-65 Hz excluded",
                 "hf_excess_over_floor": "mean PSD 70-110 / mean PSD 150-200 (descriptor: ~1 "
                                         "suggests the 70-110 band sits at a flat noise floor)",
                 "rms_*_uV": "band RMS from integrated PSD (calibration-dependent; secondary)"},
}
LINES = [60, 120, 180, 240]
HW = PARAMS["notch_exclusion_halfwidth_hz"]


# ------------------------------------------------------------------ I/O
def load_set(folder, pid):
    base = os.path.join(folder, f"sub-{pid}_task-RestingState")
    missing = [base + s for s in ("_eeg.set", "_events.tsv", "_channels.tsv") if not os.path.exists(base + s)]
    if missing:
        raise FileNotFoundError("recording files missing (Colab wipes /content on restart) -- "
                                "run setup_hbn_files.py first. Missing: " + ", ".join(os.path.basename(m) for m in missing))
    m = scipy.io.loadmat(base + "_eeg.set", squeeze_me=True, struct_as_record=False)
    if "data" not in m or isinstance(m["data"], str):
        raise RuntimeError("unsupported .set layout (expected top-level numeric 'data')")
    X = np.asarray(m["data"], dtype=np.float32)
    names = pd.read_csv(base + "_channels.tsv", sep="\t")["name"].tolist()
    if len(names) != X.shape[0]:
        raise RuntimeError(f"{len(names)} channel names vs {X.shape[0]} data rows")
    return X, float(m["srate"]), names, base + "_events.tsv"


def ec_blocks(path, sf, n):
    ev = pd.read_csv(path, sep="\t")
    ev = ev[ev["value"].isin(["instructed_toCloseEyes", "instructed_toOpenEyes"])]
    ev = ev.sort_values("onset").reset_index(drop=True)
    out = []
    for i, r in ev.iterrows():
        if r["value"] == "instructed_toCloseEyes":
            nx = ev[(ev.index > i) & (ev["value"] == "instructed_toOpenEyes")]
            end = nx["onset"].iloc[0] if len(nx) else n / sf
            s, e = int((r["onset"] + 2) * sf), int((end - 2) * sf)
            if e > s:
                out.append((s, e))
    return out


# ------------------------------------------------------------------ signal helpers
def longest_identical_run(x):
    """Length (samples) and value of the longest run of identical consecutive samples."""
    same = np.diff(x) == 0
    if not same.any():
        return 1, float(x[0])
    d = np.diff(np.concatenate([[0], same.astype(np.int8), [0]]))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    k = int(np.argmax(en - st))
    return int(en[k] - st[k] + 1), float(x[st[k]])


def preprocess_channel(x64, sf, sos_hp):
    y = sosfiltfilt(sos_hp, x64)
    for L in LINES:
        if L < sf / 2:
            b, a = iirnotch(L, 30, sf)
            y = filtfilt(b, a, y)
    return y


def _mask(f, lo, hi):
    m = (f >= lo) & (f <= hi)
    for L in LINES:
        m &= ~((f >= L - HW) & (f <= L + HW))
    return m


def spectral_features(f, P):
    def integ(lo, hi):
        m = _mask(f, lo, hi)
        return np.trapezoid(P[m], f[m]) if hasattr(np, "trapezoid") else np.trapz(P[m], f[m])
    out = {"rel_30_45": integ(30, 45) / integ(1, 45),
           "rel_70_110": integ(70, 110) / integ(1, 110),
           "hf_frac_30_110": integ(30, 110) / integ(1, 110)}
    m = (f >= 30) & (f <= 90) & ~((f >= 55) & (f <= 65)) & (P > 0)
    out["slope_30_90"] = float(np.polyfit(np.log10(f[m]), np.log10(P[m]), 1)[0]) if m.sum() > 5 else np.nan
    a, b = _mask(f, 70, 110), _mask(f, 150, 200)
    out["hf_excess_over_floor"] = float(P[a].mean() / P[b].mean()) if b.any() and P[b].mean() > 0 else np.nan
    out["rms_1_110_uV"] = float(np.sqrt(integ(1, 110)))
    out["rms_30_45_uV"] = float(np.sqrt(integ(30, 45)))
    out["rms_1_40_uV"] = float(np.sqrt(integ(1, 40)))
    return {k: float(v) for k, v in out.items()}


def fit_aperiodic(f, P, lo, hi=40):
    try:
        fm = FOOOF(peak_width_limits=[1, 8], max_n_peaks=6, aperiodic_mode="fixed", verbose=False)
        fm.fit(f, P, [lo, hi])
        return dict(exponent=float(fm.get_params("aperiodic_params", "exponent")),
                    offset=float(fm.get_params("aperiodic_params", "offset")),
                    r2=float(fm.get_params("r_squared")), error=float(fm.get_params("error")),
                    n_peaks=int(fm.n_peaks_))
    except Exception:
        return dict(exponent=np.nan, offset=np.nan, r2=np.nan, error=np.nan, n_peaks=np.nan)


def fit_aperiodic_knee(f, P, lo, hi=40):
    """Knee-mode aperiodic fit. The optimiser probes negative knees during fitting, which
    raises harmless log10 RuntimeWarnings; those are silenced here only."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            fm = FOOOF(peak_width_limits=[1, 8], max_n_peaks=6, aperiodic_mode="knee", verbose=False)
            fm.fit(f, P, [lo, hi])
        return dict(exponent=float(fm.get_params("aperiodic_params", "exponent")),
                    knee=float(fm.get_params("aperiodic_params", "knee")),
                    r2=float(fm.get_params("r_squared")))
    except Exception:
        return dict(exponent=np.nan, knee=np.nan, r2=np.nan)


def robust_z(v, idx):
    lv = np.log(np.maximum(np.asarray(v, float), 1e-30))
    med = np.nanmedian(lv[idx]); mad = np.nanmedian(np.abs(lv[idx] - med)) * 1.4826 + 1e-12
    return (lv - med) / mad


# ------------------------------------------------------------------ one participant
def process(pid, folder):
    X, sf, names, evp = load_set(folder, pid)
    n = X.shape[1]
    blocks = ec_blocks(evp, sf, n)
    if not blocks:
        raise RuntimeError("no eyes-closed blocks")
    ec = np.zeros(n, bool)
    for s, e in blocks:
        ec[s:e] = True
    sos_hp = butter(4, 1.0, btype="highpass", fs=sf, output="sos")
    sos_bs = butter(4, [7, 14], btype="bandstop", fs=sf, output="sos")
    C = len(names)
    qc = pd.DataFrame({"channel": names})
    qc["location"] = ["reference" if c == REF_NAME else ("peripheral" if c in PERIPH else "scalp")
                      for c in names]
    for g, chs in GROUPS.items():
        qc.loc[qc.channel.isin(chs), "group"] = g
    qc.loc[qc.location == "scalp", "group"] = "scalp"
    qc.loc[qc.location == "reference", "group"] = "reference"

    Y = np.zeros((C, n))
    n_nonfinite, run_s, run_val, sd_ec = (np.zeros(C), np.zeros(C), np.full(C, np.nan), np.full(C, np.nan))
    sd_ec_na = np.full(C, np.nan)
    for i in range(C):
        x = X[i]
        fin = np.isfinite(x)
        n_nonfinite[i] = (~fin).sum()
        L, v = longest_identical_run(x)
        run_s[i], run_val[i] = L / sf, v
        if fin.all():
            Y[i] = preprocess_channel(x.astype(np.float64), sf, sos_hp)
            sd_ec[i] = Y[i, ec].std()
            sd_ec_na[i] = sosfiltfilt(sos_bs, Y[i])[ec].std()          # alpha removed
    qc["n_nonfinite"], qc["longest_identical_run_s"] = n_nonfinite.astype(int), run_s
    qc["run_value"], qc["run_value_is_known_code"] = run_val, np.isclose(run_val, KNOWN_STUCK, atol=0.05)
    qc["sd_notched_hp_ec_uV"] = sd_ec
    qc["sd_notched_hp_ec_noalpha_uV"] = sd_ec_na

    def reasons(r):
        if r.location == "reference":
            return "reference_channel"
        w = []
        if r.n_nonfinite > 0: w.append("nonfinite")
        if r.longest_identical_run_s >= 1.0: w.append("stuck_run_ge_1s")
        if np.isfinite(r.sd_notched_hp_ec_uV) and r.sd_notched_hp_ec_uV < 1e-3: w.append("flat")
        return ";".join(w)
    qc["invalid_reasons"] = qc.apply(reasons, axis=1)
    qc["is_invalid"] = (qc.invalid_reasons != "") & (qc.location != "reference")

    valid_scalp = qc.index[(qc.location == "scalp") & ~qc.is_invalid].tolist()
    qc["robust_z_sd"] = np.nan
    qc.loc[valid_scalp, "robust_z_sd"] = robust_z(sd_ec_na, valid_scalp)[valid_scalp]
    qc["is_high_variance"] = (qc.location == "scalp") & (qc.robust_z_sd > 5.0)
    ref_idx = [i for i in valid_scalp if not qc.at[i, "is_high_variance"]] + [names.index(REF_NAME)]
    qc["in_reference"] = qc.index.isin(ref_idx)

    Y = Y - Y[ref_idx].mean(axis=0, keepdims=True)                  # re-reference

    # PSDs: per block, then mean
    Pb = []
    for s, e in blocks:
        f, P = welch(Y[:, s:e], fs=sf, nperseg=int(2 * sf), noverlap=int(sf), window="hann", axis=-1)
        Pb.append(P)
    Pm = np.mean(Pb, axis=0)

    # HF-elevation flag (recorded only)
    hfp = np.array([Pm[i, _mask(f, 65, 100)].sum() for i in range(C)])
    qc["robust_z_hf_65_100"] = np.nan
    qc.loc[valid_scalp, "robust_z_hf_65_100"] = robust_z(hfp, valid_scalp)[valid_scalp]
    qc["hf_elevated_recorded"] = qc.robust_z_hf_65_100 > 5.0

    # per-channel features
    rows = []
    for i, c in enumerate(names):
        r = {"participant_id": pid, **qc.iloc[i].to_dict()}
        if c != REF_NAME and not qc.at[i, "is_invalid"]:
            r.update(spectral_features(f, Pm[i]))
            a2 = fit_aperiodic(f, Pm[i], 2); a1 = fit_aperiodic(f, Pm[i], 1)
            r.update({f"ap2_{k}": v for k, v in a2.items()}); r.update({f"ap1_{k}": v for k, v in a1.items()})
            r.update({f"ap1k_{k}": v for k, v in fit_aperiodic_knee(f, Pm[i], 1).items()})
        rows.append(r)
    CH = pd.DataFrame(rows)

    # summary sets
    # "scalp" = reference-eligible channels (primary outcome, standard pipeline).
    # Contamination summaries use ALL valid scalp channels, INCLUDING high-variance ones:
    # the high-variance rule protects the reference, but must not remove contaminated
    # channels from the measurement of contamination itself.
    base_sc = CH.location.eq("scalp") & CH.in_reference & ~CH.is_invalid
    allv = CH.location.eq("scalp") & ~CH.is_invalid
    use = {"scalp": base_sc, "scalp_allvalid": allv, "frontal": allv & CH.channel.isin(FRONTAL),
           "scalp_nonfrontal": allv & ~CH.channel.isin(FRONTAL)}
    for g in GROUPS:
        use[g] = CH.group.eq(g) & ~CH.is_invalid
    rmin = PARAMS["aperiodic_fit_quality_min_r2"]
    row = dict(participant_id=pid, sfreq=sf, n_samples=n, n_ec_blocks=len(blocks),
               ec_seconds=float(ec.sum() / sf),
               n_invalid=int(qc.is_invalid.sum()), n_high_variance=int(qc.is_high_variance.sum()),
               n_reference_channels=len(ref_idx),
               invalid_channels=";".join(qc.channel[qc.is_invalid]),
               high_variance_channels=";".join(qc.channel[qc.is_high_variance]),
               hf_elevated_channels=";".join(qc.channel[qc.hf_elevated_recorded.fillna(False)]))
    feats = ["rel_30_45", "rel_70_110", "hf_frac_30_110", "slope_30_90", "hf_excess_over_floor",
             "rms_1_110_uV", "rms_30_45_uV", "rms_1_40_uV"]
    for g, m in use.items():
        row[f"{g}__n_channels"] = int(m.sum())
        for ft in feats:
            row[f"{g}__{ft}"] = float(CH.loc[m, ft].median()) if m.any() else np.nan
    def med(mask, col):
        return float(CH.loc[mask, col].median()) if mask.any() else np.nan
    for g in ("scalp", "scalp_allvalid", "frontal", "scalp_nonfrontal"):
        row[f"{g}__n_high_variance"] = int((use[g] & CH.is_high_variance).sum())
        row[f"{g}__n_hf_elevated"] = int((use[g] & CH.hf_elevated_recorded.fillna(False).astype(bool)).sum())
        for tag in ("ap2", "ap1"):
            ok = use[g] & (CH[f"{tag}_r2"] >= rmin)
            nall = int(use[g].sum())
            row[f"{g}__{tag}_n_fit_ok"] = int(ok.sum())
            row[f"{g}__{tag}_n_fit_fail"] = nall - int(ok.sum())
            row[f"{g}__{tag}_frac_fit_fail"] = (nall - int(ok.sum())) / nall if nall else np.nan
            row[f"{g}__{tag}_exponent"] = med(ok, f"{tag}_exponent")          # R^2-filtered
            row[f"{g}__{tag}_exponent_nofilter"] = med(use[g], f"{tag}_exponent")
            # mean is deliberately NOT robust: a median ignores a contaminated minority
            row[f"{g}__{tag}_exponent_nofilter_mean"] = (float(CH.loc[use[g], f"{tag}_exponent"].mean())
                                                          if use[g].any() else np.nan)
            row[f"{g}__{tag}_offset"] = med(ok, f"{tag}_offset")
            row[f"{g}__{tag}_median_r2"] = med(use[g], f"{tag}_r2")
    okk = use["scalp"] & (CH["ap1k_r2"] >= rmin)
    row["scalp__ap1k_n_fit_ok"] = int(okk.sum())
    row["scalp__ap1k_exponent"] = med(okk, "ap1k_exponent")
    row["scalp__ap1k_knee"] = med(okk, "ap1k_knee")
    row["scalp__ap1k_median_r2"] = med(use["scalp"], "ap1k_r2")
    for g in GROUPS:
        row[f"gradient_{g}_vs_scalp__rel_70_110"] = row[f"{g}__rel_70_110"] / row["scalp__rel_70_110"]
    for ft in ("rel_70_110", "hf_frac_30_110", "rel_30_45"):
        row[f"gradient_frontal_vs_nonfrontal__{ft}"] = row[f"frontal__{ft}"] / row[f"scalp_nonfrontal__{ft}"]
    # Median gradients are blind to muscle confined to a MINORITY of channels. This metric is
    # not: the fraction of a group's channels exceeding the 90th percentile of non-frontal
    # scalp. Its null expectation (no localized excess) is ~0.10 by construction.
    for ft in ("rel_70_110", "hf_frac_30_110"):
        ref = CH.loc[use["scalp_nonfrontal"], ft].dropna()
        thr = float(np.percentile(ref, 90)) if len(ref) else np.nan
        row[f"scalp_nonfrontal__{ft}_p90"] = thr
        for g in ("infraorbital", "face_jaw", "neck", "frontal"):
            v = CH.loc[use[g], ft].dropna()
            row[f"{g}__{ft}_frac_above_nonfrontal_p90"] = float((v > thr).mean()) if len(v) else np.nan
            row[f"{g}__{ft}_mean"] = float(v.mean()) if len(v) else np.nan

    # block-level stability: primary exponent + primary spectral features
    sc_idx = CH.index[use["scalp"]].tolist()
    brows = []
    for b, P in enumerate(Pb):
        ex = [fit_aperiodic(f, P[i], 2) for i in sc_idx]
        ex = [d["exponent"] for d in ex if d["r2"] >= rmin]
        br = dict(participant_id=pid, block=b + 1, scalp_ap2_exponent=float(np.median(ex)) if ex else np.nan)
        for g, m in use.items():
            idx = CH.index[m].tolist()
            vals = pd.DataFrame([spectral_features(f, P[i]) for i in idx])
            for ft in ("rel_70_110", "hf_frac_30_110", "slope_30_90"):
                br[f"{g}__{ft}"] = float(vals[ft].median()) if len(vals) else np.nan
        brows.append(br)
    BL = pd.DataFrame(brows)
    for col in [c for c in BL.columns if c not in ("participant_id", "block")]:
        v = BL[col].dropna()
        row[f"{col}__block_sd"] = float(v.std(ddof=1)) if len(v) > 1 else np.nan
        if col.endswith("exponent"):                   # CV only where the quantity is strictly positive
            row[f"{col}__block_cv"] = float(v.std(ddof=1) / v.mean()) if len(v) > 1 and v.mean() > 0 else np.nan
    return row, CH, BL


def run(participants, out):
    os.makedirs(out, exist_ok=True)
    P, Cs, Bs = [], [], []
    for pid, folder in participants.items():
        try:
            row, CH, BL = process(pid, folder)
            row.update(excluded=False, exclusion_reason="")
            Cs.append(CH); Bs.append(BL)
            print(f"  OK   {pid}: invalid={row['n_invalid']} high-var={row['n_high_variance']} "
                  f"ref={row['n_reference_channels']} | scalp exponent (2-40) = "
                  f"{row['scalp__ap2_exponent']:.3f} from {row['scalp__ap2_n_fit_ok']} channels")
        except Exception as e:
            row = dict(participant_id=pid, excluded=True, exclusion_reason=str(e))
            print(f"  FAIL {pid}: {e}")
        P.append(row)
    pd.DataFrame(P).to_csv(f"{out}/v31_participants.csv", index=False)
    if Cs:
        pd.concat(Cs).to_csv(f"{out}/v31_channels.csv", index=False)
        pd.concat(Bs).to_csv(f"{out}/v31_blocks.csv", index=False)
    prov = dict(params=PARAMS, groups=GROUPS, python=platform.python_version(),
                numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__,
                fooof=fooof.__version__, participants=participants)
    with open(f"{out}/v31_provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2, default=str)
    print(f"\nWrote {out}/v31_participants.csv, v31_channels.csv, v31_blocks.csv, v31_provenance.json")
    return pd.DataFrame(P)


if __name__ == "__main__":          # pasted as a Colab cell -> runs the 3-participant pilot
    summary = run(PARTICIPANTS, OUT)
