# =============================================================================
# HBN R3 -- PRE-REGISTERED ANALYSIS for pipeline v3.1      (Colab cell)
# Implements "Analysis Plan -- Peripheral High-Frequency Activity and the Scalp
# Aperiodic Exponent (HBN Release 3)", fixed 20 October 2026, with deviation D1
# (criterion 2 redefined before any cohort data was inspected) and deviation D2
# (coding error in the fallback attention analysis corrected). Every number below
# comes from that plan; change NOTHING here without a dated entry in its
# deviations log.
#
# Flow: load per-participant results -> sample flow + exclusions -> Q1 gatekeeper
#       -> (pass) Q2 primary + Q3 secondary   |   (fail) fallback
#       -> sensitivity analyses -> report files in OUT/analysis/
# Refuses to run on an incomplete cohort unless ALLOW_PARTIAL = True (testing only).
# =============================================================================
import os, json, glob, hashlib, warnings
import numpy as np, pandas as pd
from scipy import stats
try:
    import statsmodels.api as sm, patsy
except ImportError as e:
    raise ImportError("statsmodels/patsy missing (both ship with Colab)") from e

OUT = "/content/drive/MyDrive/hbn/R3_v31"
INCLUDE_STATUS = ("available", "caution")
ALLOW_PARTIAL = False            # True ONLY for testing; the plan analyses the complete cohort

# ---- plan constants (Analysis Plan, 2026-10-20) --------------------------------
N_BOOT, SEED = 2000, 20261020
REL_MIN = 0.70                   # Q1 criterion 1: split-half reliability (Spearman-Brown)
MIN_SCALP_CH, MIN_EC_S = 80, 120 # exclusions
HALF_A, HALF_B = (1, 3, 5), (2, 4)
SPLINE = "cr(age, df=3, constraints='center')"   # centred: an uncentred basis is collinear
                                                 # with the intercept (verified)


# ------------------------------------------------------------------ load + guard
def load():
    expected = pd.read_csv(os.path.join(OUT, "participants.tsv"), sep="\t")
    expected["pid"] = expected["participant_id"].str.replace("^sub-", "", regex=True)
    expected = expected[expected["RestingState"].isin(INCLUDE_STATUS)]
    base = os.path.join(OUT, "per_participant")
    P, B, failed, h = [], [], [], hashlib.sha256()
    for pid in sorted(expected["pid"]):
        d = os.path.join(base, pid)
        if os.path.exists(os.path.join(d, "participant.csv")):
            for f, store in (("participant.csv", P), ("blocks.csv", B)):
                with open(os.path.join(d, f), "rb") as fh:
                    h.update(fh.read())
                store.append(pd.read_csv(os.path.join(d, f)))
        elif os.path.exists(os.path.join(d, "failed.json")):
            with open(os.path.join(d, "failed.json")) as fh:
                failed.append(json.load(fh))
    pending = len(expected) - len(P) - len(failed)
    if pending and not ALLOW_PARTIAL:
        raise RuntimeError(f"{pending} participants not yet processed. Finish the cohort run first "
                           "(the plan analyses the complete cohort).")
    return expected, pd.concat(P, ignore_index=True), pd.concat(B, ignore_index=True), failed, pending, h.hexdigest()


def log10p(x):
    x = np.asarray(x, float)
    return np.where(x > 0, np.log10(np.where(x > 0, x, 1.0)), np.nan)


def z(x):
    x = np.asarray(x, float)
    return (x - np.nanmean(x)) / np.nanstd(x, ddof=1)


# ------------------------------------------------------------------ statistics
def boot_ci(stat, n, rng):
    """Percentile CI over N_BOOT resamples of children (rows). stat(idx) -> float or array."""
    vals, fails = [], 0
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        try:
            v = stat(idx)
            if np.all(np.isfinite(v)):
                vals.append(v)
            else:
                fails += 1
        except Exception:
            fails += 1
    vals = np.array(vals)
    return np.percentile(vals, 2.5, axis=0), np.percentile(vals, 97.5, axis=0), fails


def partial_rank_corr(x, y, covs):
    """Partial Spearman: rank everything, residualise ranks of x and y on ranks of covs."""
    R = lambda v: stats.rankdata(v)
    C = np.column_stack([np.ones(len(x))] + [R(c) for c in covs])
    rx = R(x) - C @ np.linalg.lstsq(C, R(x), rcond=None)[0]
    ry = R(y) - C @ np.linalg.lstsq(C, R(y), rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def huber_fit(y, X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit().params


def design(df, terms):
    """Design matrix built ONCE on the analysis sample (spline knots fixed); bootstrap
    resamples its rows."""
    # isolated namespace: names left in a shared notebook (e.g. a table called C) must not
    # shadow patsy built-ins such as C() and cr()
    return patsy.dmatrix("1 + " + " + ".join(terms), df, return_type="dataframe",
                         eval_env=patsy.EvalEnvironment([{}]))


def robust_beta(df, outcome, predictor, extra=(), rng=None):
    d = df[[outcome, predictor, "age", "sex", *extra]].dropna().copy()
    d["_y"], d["_x"] = z(d[outcome]), z(d[predictor])
    X = design(d, ["_x", SPLINE, "C(sex)", *extra])
    y, Xv, j = d["_y"].values, X.values, list(X.columns).index("_x")
    beta = float(huber_fit(y, Xv)[j])
    lo, hi, fails = boot_ci(lambda i: huber_fit(y[i], Xv[i])[j], len(y), rng)
    return dict(n=len(y), beta=beta, ci_low=float(lo), ci_high=float(hi),
                excludes_zero=bool(lo > 0 or hi < 0), boot_failures=fails)


# ------------------------------------------------------------------ analysis
def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(os.path.join(OUT, "analysis"), exist_ok=True)
    expected, P, B, failed, pending, digest = load()
    R = {"plan_date": "2026-10-20", "input_sha256": digest, "seed": SEED, "n_boot": N_BOOT}

    # ---- sample flow + exclusions
    flow = {"listed (status included)": len(expected), "processed": len(P),
            "failed processing": len(failed), "not yet processed": pending}
    m_ch = P["scalp__n_channels"] >= MIN_SCALP_CH
    m_ec = P["ec_seconds"] >= MIN_EC_S
    flow[f"excluded: < {MIN_SCALP_CH} scalp channels"] = int((~m_ch).sum())
    flow[f"excluded: < {MIN_EC_S} s eyes closed (and >= {MIN_SCALP_CH} ch)"] = int((m_ch & ~m_ec).sum())
    S = P[m_ch & m_ec].copy()
    flow["analysed"] = len(S)
    R["flow"] = flow

    # ---- derived variables (plan: Measures)
    S["PHF"] = (log10p(S["face_jaw__rel_70_110"]) + log10p(S["neck__rel_70_110"])) / 2
    S["PHF_neck"] = log10p(S["neck__rel_70_110"])
    S["PHF_jaw"] = log10p(S["face_jaw__rel_70_110"])
    S["PHF_3045"] = (log10p(S["face_jaw__rel_30_45"]) + log10p(S["neck__rel_30_45"])) / 2
    S["log_frontal"] = log10p(S["frontal__rel_70_110"])
    S["log_nonfrontal"] = log10p(S["scalp_nonfrontal__rel_70_110"])
    S["frontal_minus_nonfrontal"] = S["log_frontal"] - S["log_nonfrontal"]   # deviation D1
    S["log_floor"] = log10p(S["scalp__hf_excess_over_floor"])
    S["D"] = S["scalp_allvalid__ap2_exponent_nofilter_mean"] - S["scalp__ap2_exponent"]
    S["D_ap1"] = S["scalp_allvalid__ap1_exponent_nofilter_mean"] - S["scalp__ap1_exponent"]

    # ---- Q1 gatekeeper
    Bs = B[B["participant_id"].isin(S["participant_id"])].copy()
    Bs["PHF_b"] = (log10p(Bs["face_jaw__rel_70_110"]) + log10p(Bs["neck__rel_70_110"])) / 2
    ha = Bs[Bs.block.isin(HALF_A)].groupby("participant_id")["PHF_b"].mean()
    hb = Bs[Bs.block.isin(HALF_B)].groupby("participant_id")["PHF_b"].mean()
    H = pd.concat([ha.rename("a"), hb.rename("b")], axis=1).dropna()
    sb = lambda r: 2 * r / (1 + r)
    r_half = float(np.corrcoef(H["a"], H["b"])[0, 1])
    lo, hi, _ = boot_ci(lambda i: sb(np.corrcoef(H["a"].values[i], H["b"].values[i])[0, 1]), len(H), rng)
    crit1 = dict(n=len(H), r_half=r_half, spearman_brown=sb(r_half), ci_low=float(lo), ci_high=float(hi),
                 threshold=REL_MIN, passes=bool(sb(r_half) >= REL_MIN))
    # Deviation D1: Spearman rho(PHF, frontal-minus-nonfrontal log HF). A scalp-wide floor
    # cancels in the within-person difference; the original partial correlation did not
    # (controlling a NOISY floor proxy left residual confounding: it passed a pure floor).
    g = S[["PHF", "frontal_minus_nonfrontal"]].dropna().values
    rho = float(stats.spearmanr(g[:, 0], g[:, 1])[0])
    lo, hi, _ = boot_ci(lambda i: stats.spearmanr(g[i, 0], g[i, 1])[0], len(g), rng)
    crit2 = dict(n=len(g), spearman_rho_PHF_vs_frontal_minus_nonfrontal=rho, ci_low=float(lo),
                 ci_high=float(hi), passes=bool(lo > 0),
                 limitation="cannot exclude a noise floor that loads more strongly on frontal channels")
    gate = crit1["passes"] and crit2["passes"]
    R["Q1"] = dict(criterion1_reliability=crit1, criterion2_spatial_specificity=crit2, passes=gate)

    sens = []
    if gate:
        # ---- Q2 primary
        q2 = robust_beta(S, "scalp__ap2_exponent", "PHF", rng=rng)
        d = S[["scalp__ap2_exponent", "PHF", "age", "sex"]].dropna()
        q2["spearman_rho"] = float(stats.spearmanr(d["PHF"], d["scalp__ap2_exponent"])[0])
        q2["partial_spearman_age_sex"] = partial_rank_corr(
            d["PHF"].values, d["scalp__ap2_exponent"].values,
            [d["age"].values, (d["sex"] == d["sex"].iloc[0]).astype(float).values])
        q2["significant"] = q2["excludes_zero"]
        R["Q2_primary"] = q2

        # ---- Q3 secondary (estimation only)
        d = S[["scalp__ap2_exponent", "PHF", "attention", "age", "sex"]].dropna().copy()
        d["_y"], d["_a"], d["_p"] = z(d["scalp__ap2_exponent"]), z(d["attention"]), z(d["PHF"])
        XA, XB = design(d, ["_a", SPLINE, "C(sex)"]), design(d, ["_a", "_p", SPLINE, "C(sex)"])
        y, A, Bm = d["_y"].values, XA.values, XB.values
        ja, jb = list(XA.columns).index("_a"), list(XB.columns).index("_a")
        bA, bB = float(huber_fit(y, A)[ja]), float(huber_fit(y, Bm)[jb])
        lo, hi, _ = boot_ci(lambda i: huber_fit(y[i], A[i])[ja] - huber_fit(y[i], Bm[i])[jb], len(y), rng)
        R["Q3_secondary"] = dict(n=len(y), beta_attention_A=bA, beta_attention_B=bB,
                                 delta=bA - bB, delta_ci_low=float(lo), delta_ci_high=float(hi))

        # ---- sensitivity (descriptive)
        specs = [("1 exclude caution", S[S["curator_resting_status"] != "caution"], "scalp__ap2_exponent", "PHF", ()),
                 ("2 exponent 1-40 Hz", S, "scalp__ap1_exponent", "PHF", ()),
                 ("3 exponent knee mode", S, "scalp__ap1k_exponent", "PHF", ()),
                 ("4 uncensored exponent", S, "scalp_allvalid__ap2_exponent_nofilter_mean", "PHF", ()),
                 ("5a PHF neck only", S, "scalp__ap2_exponent", "PHF_neck", ()),
                 ("5b PHF face/jaw only", S, "scalp__ap2_exponent", "PHF_jaw", ()),
                 ("6 PHF 30-45 Hz", S, "scalp__ap2_exponent", "PHF_3045", ()),
                 ("7 + noise-floor covariate", S, "scalp__ap2_exponent", "PHF", ("log_floor",))]
        for name, df, yv, xv, ex in specs:
            sens.append(dict(analysis=name, **robust_beta(df, yv, xv, ex, rng=rng)))
    else:
        # ---- fallback
        Dv = S["D"].dropna().values
        w = stats.wilcoxon(Dv)
        lo, hi, _ = boot_ci(lambda i: np.median(Dv[i]), len(Dv), rng)
        fb = dict(n=len(Dv), median_D=float(np.median(Dv)), mean_D=float(np.mean(Dv)),
                  wilcoxon_p=float(w.pvalue), ci_low=float(lo), ci_high=float(hi),
                  significant=bool(w.pvalue < 0.05))
        d = S[["D", "log_frontal"]].dropna().values
        rho = float(stats.spearmanr(d[:, 0], d[:, 1])[0])
        lo, hi, _ = boot_ci(lambda i: stats.spearmanr(d[i, 0], d[i, 1])[0], len(d), rng)
        fb["secondary_rho_D_frontalHF"] = dict(rho=rho, ci_low=float(lo), ci_high=float(hi))
        d = S[["scalp__ap2_exponent", "scalp_allvalid__ap2_exponent_nofilter_mean", "attention", "age", "sex"]].dropna().copy()
        d["_a"] = z(d["attention"]); Xdf = design(d, ["_a", SPLINE, "C(sex)"]); X = Xdf.values
        # Deviation D2 (coding error corrected): this was hard-coded j = 1, but patsy places
        # categorical terms first, so column 1 is C(sex)[T.M]. Look the column up by name.
        j = list(Xdf.columns).index("_a")
        ys, yu = z(d["scalp__ap2_exponent"]), z(d["scalp_allvalid__ap2_exponent_nofilter_mean"])
        bs_, bu_ = float(huber_fit(ys, X)[j]), float(huber_fit(yu, X)[j])
        lo, hi, _ = boot_ci(lambda i: huber_fit(ys[i], X[i])[j] - huber_fit(yu[i], X[i])[j], len(ys), rng)
        fb["secondary_attention"] = dict(n=len(ys), beta_standard=bs_, beta_uncensored=bu_,
                                         difference=bs_ - bu_, ci_low=float(lo), ci_high=float(hi))
        fb["descriptive"] = dict(median_fit_failures=float(S["scalp_allvalid__ap2_n_fit_fail"].median()),
                                 median_high_variance=float(S["n_high_variance"].median()))
        R["fallback"] = fb
        for name, df, col in [("1 exclude caution", S[S["curator_resting_status"] != "caution"], "D"),
                              ("2 exponent 1-40 Hz", S, "D_ap1")]:
            v = df[col].dropna().values
            lo, hi, _ = boot_ci(lambda i: np.median(v[i]), len(v), rng)
            sens.append(dict(analysis=name, n=len(v), median_D=float(np.median(v)),
                             wilcoxon_p=float(stats.wilcoxon(v).pvalue), ci_low=float(lo), ci_high=float(hi)))
        for name in ("3 knee mode", "4 uncensored exponent", "5 PHF variants", "6 PHF 30-45 Hz", "7 noise floor"):
            sens.append(dict(analysis=name, note="not applicable to the fallback branch (plan)"))
    R["sensitivity"] = sens

    # ---- write + print
    A = os.path.join(OUT, "analysis")
    with open(os.path.join(A, "results.json"), "w") as fh:
        json.dump(R, fh, indent=2, default=float)
    pd.DataFrame(sens).to_csv(os.path.join(A, "sensitivity.csv"), index=False)
    pd.DataFrame([flow]).T.rename(columns={0: "n"}).to_csv(os.path.join(A, "sample_flow.csv"))
    S.to_csv(os.path.join(A, "analysis_sample.csv"), index=False)
    print("SAMPLE FLOW"); [print(f"  {k:<48}{v}") for k, v in flow.items()]
    c1, c2 = R["Q1"]["criterion1_reliability"], R["Q1"]["criterion2_spatial_specificity"]
    print(f"\nQ1 GATEKEEPER: {'PASS' if gate else 'FAIL'}")
    print(f"  reliability (Spearman-Brown) {c1['spearman_brown']:.3f}  [{c1['ci_low']:.3f}, {c1['ci_high']:.3f}]"
          f"  threshold {REL_MIN} -> {'pass' if c1['passes'] else 'fail'}")
    print(f"  rho PHF ~ (frontal - nonfrontal) {c2['spearman_rho_PHF_vs_frontal_minus_nonfrontal']:.3f}  [{c2['ci_low']:.3f}, {c2['ci_high']:.3f}]"
          f" -> {'pass' if c2['passes'] else 'fail'}")
    if gate:
        q = R["Q2_primary"]
        print(f"\nQ2 PRIMARY  beta_PHF = {q['beta']:.3f}  95% CI [{q['ci_low']:.3f}, {q['ci_high']:.3f}]  n={q['n']}"
              f"  -> {'significant' if q['significant'] else 'not significant'}")
        q = R["Q3_secondary"]
        print(f"Q3 SECONDARY beta_attention {q['beta_attention_A']:.3f} -> {q['beta_attention_B']:.3f}"
              f"  delta {q['delta']:.3f} [{q['delta_ci_low']:.3f}, {q['delta_ci_high']:.3f}]  n={q['n']}")
    else:
        q = R["fallback"]
        print(f"\nFALLBACK PRIMARY median D = {q['median_D']:.4f}  95% CI [{q['ci_low']:.4f}, {q['ci_high']:.4f}]"
              f"  Wilcoxon p = {q['wilcoxon_p']:.3g}  n={q['n']}")
    print(f"\nSENSITIVITY (descriptive):"); print(pd.DataFrame(sens).round(3).to_string(index=False))
    print(f"\nWritten to {A}/  (input sha256 {digest[:16]}...)")
    return R


if __name__ == "__main__":
    results = main()
