# =============================================================================
# REVISION P0-1, P0-2, P0-4      (Colab cell; implements revision_P0_spec.md, approved)
#
# Read-only on all existing files. Writes ONLY to a new timestamped folder:
#   /content/drive/MyDrive/hbn/revision_P0/run_<YYYYmmdd_HHMMSS>/
# No new model fits: everything is recomputed from per-channel fits already stored.
# Estimation only (bootstrap CIs); no new p-values. Gates G1-G3 stop the run on failure.
# Run after: drive.mount('/content/drive')
# =============================================================================
import os, sys, json, glob, time, struct, hashlib, platform, warnings, itertools
import numpy as np, pandas as pd, scipy
from scipy import stats
import statsmodels, statsmodels.api as sm, patsy

DRIVE    = "/content/drive/MyDrive/hbn"
R3_DIR   = f"{DRIVE}/R3_v31"
SCAN_DIR = f"{DRIVE}/stuck_scan/per_participant"
OUT_BASE = f"{DRIVE}/revision_P0"
RUN_STAMP = None                     # None = current time; the folder must not already exist
N_BOOT, SEED, RMIN = 2000, 20261022, 0.85
TOL_G1, TOL_G2 = 1e-9, 1e-12
SCAN_RELEASES = ("R1", "R2", "R4", "R5", "R6", "R7", "R8", "R9")
PERIPH = {"E48","E119","E127","E126","E113","E49","E56","E107","E99","E63",
          "E94","E68","E43","E120","E88","E73","E81","E125","E128"}
SPLINE = "cr(age, df=3, constraints='center')"


class GateFailure(RuntimeError):
    pass


def _bool(s):
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1"))


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _boot(stat, n, rng):
    vals, fails = [], 0
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        try:
            v = stat(i)
            if np.all(np.isfinite(v)): vals.append(v)
            else: fails += 1
        except Exception:
            fails += 1
    vals = np.array(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), fails


def _design(df, formula):
    # isolated namespace: names left in the notebook cannot shadow C() or cr()
    return patsy.dmatrix(formula, df, return_type="dataframe", eval_env=patsy.EvalEnvironment([{}]))


def _huber(y, X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit().params


def _icc(x1, x2):
    """Shrout-Fleiss two-way ICCs for n subjects x k=2 measurements."""
    Y = np.column_stack([x1, x2]); n, k = Y.shape; gm = Y.mean()
    ssr = k * ((Y.mean(1) - gm) ** 2).sum(); ssc = n * ((Y.mean(0) - gm) ** 2).sum()
    sse = ((Y - gm) ** 2).sum() - ssr - ssc
    msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
    icc31 = (msr - mse) / (msr + (k - 1) * mse)
    icc21 = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    return icc31, icc21


def _desc(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return dict(n=int(len(v)), median=float(np.median(v)), mean=float(np.mean(v)), sd=float(np.std(v, ddof=1)),
                q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)),
                p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)))


# ----------------------------------------------------------------------------- P0-1
def per_child_factorial(CH):
    """Eight summaries: stat {median, mean} x set {R, V} x filter {F, none}, 2-40 Hz fit."""
    out = []
    for pid, g in CH.groupby("participant_id", sort=True):
        V = g.location.eq("scalp") & ~g.is_invalid
        R = V & g.in_reference
        F = g.ap2_r2 >= RMIN                                  # NaN R^2 fails the filter
        ex = g.ap2_exponent
        row = dict(participant_id=pid, n_V=int(V.sum()), n_R=int(R.sum()),
                   n_hv_V=int((V & g.is_high_variance).sum()), n_fail_V=int((V & ~F).sum()),
                   n_nanfit_V=int((V & ex.isna()).sum()))
        for st in ("median", "mean"):
            for sname, smask in (("R", R), ("V", V)):
                for fname, fmask in (("F", F), ("none", pd.Series(True, index=g.index))):
                    vals = ex[smask & fmask].dropna()
                    row[f"{st}_{sname}_{fname}"] = (float(vals.median()) if st == "median" else float(vals.mean())) \
                        if len(vals) else np.nan
        out.append(row)
    return pd.DataFrame(out)


def gate_g1(T, S):
    checks = {"median_R_F": "scalp__ap2_exponent", "median_R_none": "scalp__ap2_exponent_nofilter",
              "median_V_none": "scalp_allvalid__ap2_exponent_nofilter",
              "mean_V_none": "scalp_allvalid__ap2_exponent_nofilter_mean"}
    M = T.merge(S[["participant_id", "D"] + list(checks.values())], on="participant_id", how="inner")
    if len(M) != len(S):
        raise GateFailure(f"G1: {len(S) - len(M)} sample children have no per-channel rows")
    res = {}
    pairs = list(checks.items()) + [("_D_recomputed", "D")]
    M["_D_recomputed"] = M["mean_V_none"] - M["median_R_F"]
    for mine, saved in pairs:
        a, b = M[mine].values, M[saved].values
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.where(both_nan, 0.0, np.abs(a - b))
        bad = int(np.sum(~both_nan & ~(diff <= TOL_G1)))
        res[f"{mine} vs {saved}"] = dict(max_abs_diff=float(np.nanmax(diff)) if len(diff) else 0.0, n_mismatch=bad)
    if any(v["n_mismatch"] for v in res.values()):
        raise GateFailure("G1 reproduction check failed: " + json.dumps(res))
    return res


def orderings(row):
    """Contribution of each factor under each of the 6 orders, from start cell S to end cell U."""
    start = dict(stat="median", set="R", filt="F")
    end = dict(stat="mean", set="V", filt="none")
    cell = lambda c: row[f"{c['stat']}_{c['set']}_{c['filt']}"]
    res = {}
    for order in itertools.permutations(("filt", "set", "stat")):
        cur = dict(start); prev = cell(cur); contrib = {}
        for f in order:
            cur[f] = end[f]; now = cell(cur); contrib[f] = now - prev; prev = now
        res[">".join(order)] = contrib
    return res


# ----------------------------------------------------------------------------- main
def main():
    if OUT_BASE.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    stamp = RUN_STAMP or time.strftime("%Y%m%d_%H%M%S")
    RUN = os.path.join(OUT_BASE, f"run_{stamp}")
    os.makedirs(RUN, exist_ok=False)                      # never overwrite an earlier run
    rng = np.random.default_rng(SEED)

    paths = {"analysis_sample": f"{R3_DIR}/analysis/analysis_sample.csv",
             "channels": f"{R3_DIR}/R3_v31_channels.csv", "participants": f"{R3_DIR}/R3_v31_participants.csv"}
    scan_files = sorted(f for r in SCAN_RELEASES for f in glob.glob(f"{SCAN_DIR}/{r}_*.csv"))
    log = dict(spec="revision_P0_spec.md (approved)", seed=SEED, n_boot=N_BOOT, r2_min=RMIN,
               inputs={k: dict(path=p, sha256=_sha(p)) for k, p in paths.items()},
               scan_files=dict(n=len(scan_files), sha256_of_concatenated=hashlib.sha256(
                   b"".join(open(f, "rb").read() for f in scan_files)).hexdigest()),
               versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                             pandas=pd.__version__, statsmodels=statsmodels.__version__, patsy=patsy.__version__),
               gates={})
    write_log = lambda: json.dump(log, open(os.path.join(RUN, "run_log.json"), "w"), indent=2, default=str)

    S = pd.read_csv(paths["analysis_sample"])
    CH = pd.read_csv(paths["channels"])
    for c in ("is_invalid", "is_high_variance", "in_reference"):
        CH[c] = _bool(CH[c])
    CH_s = CH[CH.participant_id.isin(S.participant_id)].copy()

    try:
        # ------------------------------------------------------------ P0-1
        T = per_child_factorial(CH_s)
        log["gates"]["G1"] = gate_g1(T, S)
        T = T.merge(S[["participant_id", "D"]], on="participant_id")
        T["S"], T["B"], T["A"], T["U"] = T.median_R_F, T.median_R_none, T.median_V_none, T.mean_V_none
        T["C"], T["H"], T["G"] = T.B - T.S, T.A - T.B, T.U - T.A
        g2 = np.abs((T.C + T.H + T.G) - T.D).max()
        log["gates"]["G2"] = dict(max_abs_diff=float(g2))
        if not g2 <= TOL_G2: raise GateFailure(f"G2: C+H+G differs from D by {g2}")
        T["subgroup"] = np.where((T.n_hv_V == 0) & (T.n_fail_V == 0), "nothing_excluded", "at_least_one_excluded")
        ne = T[T.subgroup == "nothing_excluded"]
        g3 = float(np.abs(ne.D - ne.G).max()) if len(ne) else 0.0
        log["gates"]["G3"] = dict(n_nothing_excluded=int(len(ne)), max_abs_diff=g3)
        if not g3 <= TOL_G2: raise GateFailure(f"G3: in nothing-excluded children D != G (max {g3})")

        O = [orderings(r) for _, r in T.iterrows()]
        for o in O[0]:
            for f in ("filt", "set", "stat"):
                T[f"ord_{o}_{f}"] = [x[o][f] for x in O]
        for f in ("filt", "set", "stat"):
            T[f"shapley_{f}"] = T[[f"ord_{o}_{f}" for o in O[0]]].mean(axis=1)
        T.to_csv(os.path.join(RUN, "P0-1_per_child.csv"), index=False)

        def summarise(col, label):
            v = T[col].dropna().values
            lo, hi, fl = _boot(lambda i: np.median(v[i]), len(v), rng)
            return dict(term=label, column=col, n=len(v), median=float(np.median(v)), ci_low=lo, ci_high=hi,
                        mean=float(v.mean()), sd=float(v.std(ddof=1)),
                        pct_exactly_zero=float(100 * np.mean(np.abs(v) <= TOL_G2)), boot_failures=fl)
        rows = [summarise("D", "D (total)"), summarise("C", "C = B - S (R^2 censoring, primary order)"),
                summarise("H", "H = A - B (high-variance re-inclusion, primary order)"),
                summarise("G", "G = U - A (mean vs median, primary order)")]
        names = {"filt": "R^2 filter", "set": "channel set", "stat": "statistic"}
        for f in ("filt", "set", "stat"):
            rows.append(summarise(f"shapley_{f}", f"Shapley: {names[f]}"))
        for o in O[0]:
            for f in ("filt", "set", "stat"):
                rows.append(summarise(f"ord_{o}_{f}", f"order {o}: {names[f]}"))
        pd.DataFrame(rows).to_csv(os.path.join(RUN, "P0-1_summary.csv"), index=False)
        sg = []
        for gname, gg in T.groupby("subgroup"):
            v = gg.D.dropna().values
            lo, hi, _ = _boot(lambda i: np.median(v[i]), len(v), rng)
            sg.append(dict(subgroup=gname, n=len(v), median_D=float(np.median(v)), ci_low=lo, ci_high=hi,
                           mean_D=float(v.mean())))
        pd.DataFrame(sg).to_csv(os.path.join(RUN, "P0-1_subgroups.csv"), index=False)

        # ------------------------------------------------------------ P0-2
        pc = []
        for pid, g in CH_s.groupby("participant_id", sort=True):
            R = g.location.eq("scalp") & ~g.is_invalid & g.in_reference
            m = R & (g.ap2_r2 >= RMIN) & (g.ap1_r2 >= RMIN) & g.ap2_exponent.notna() & g.ap1_exponent.notna()
            pc.append(dict(participant_id=pid, n_matched=int(m.sum()),
                           E2m=float(g.ap2_exponent[m].median()) if m.any() else np.nan,
                           E1m=float(g.ap1_exponent[m].median()) if m.any() else np.nan))
        P2 = pd.DataFrame(pc).merge(S[["participant_id", "scalp__ap2_exponent", "scalp__ap1_exponent",
                                       "age", "sex", "attention", "log_frontal"]], on="participant_id")
        P2 = P2.rename(columns={"scalp__ap2_exponent": "E2", "scalp__ap1_exponent": "E1"})
        P2["delta_matched"], P2["delta_reported"] = P2.E1m - P2.E2m, P2.E1 - P2.E2
        P2.to_csv(os.path.join(RUN, "P0-2_per_child.csv"), index=False)
        sd_ref = float(S["scalp__ap2_exponent"].std(ddof=1))
        R2 = dict(sd_units_reference=dict(column="scalp__ap2_exponent", sd=sd_ref),
                  n_empty_matched_set=int((P2.n_matched == 0).sum()), versions={})
        for ver, e1, e2, dl in (("matched_primary", "E1m", "E2m", "delta_matched"),
                                ("as_reported_secondary", "E1", "E2", "delta_reported")):
            out = {}
            d = P2[[e1, e2, dl, "age", "sex", "attention", "log_frontal"]].dropna(subset=[e1, e2])
            out["delta"] = _desc(d[dl]); out["delta_sd_units"] = {k: (v / sd_ref if k != "n" else v)
                                                                  for k, v in _desc(d[dl]).items()}
            a1, a2 = d[e1].values, d[e2].values
            rho = stats.spearmanr(a1, a2)[0]; i31, i21 = _icc(a1, a2)
            out["agreement"] = dict(
                spearman=dict(est=float(rho), ci=_boot(lambda i: stats.spearmanr(a1[i], a2[i])[0], len(a1), rng)[:2]),
                icc_3_1=dict(est=float(i31), ci=_boot(lambda i: _icc(a1[i], a2[i])[0], len(a1), rng)[:2]),
                icc_2_1=dict(est=float(i21), ci=_boot(lambda i: _icc(a1[i], a2[i])[1], len(a1), rng)[:2]), n=len(a1))
            out["delta_vs_predictors"] = {}
            for pred in ("age", "log_frontal", "attention"):
                q = d[[dl, pred]].dropna(); x, y = q[dl].values, q[pred].values
                out["delta_vs_predictors"][pred] = dict(
                    n=len(q), spearman=float(stats.spearmanr(x, y)[0]),
                    ci=_boot(lambda i: stats.spearmanr(x[i], y[i])[0], len(x), rng)[:2])
            # age slope, paired across fit ranges
            q = d[[e1, e2, "age", "sex"]].dropna().copy()
            X = _design(q, "1 + age + C(sex)"); j = list(X.columns).index("age"); Xv = X.values
            y1, y2 = q[e1].values, q[e2].values
            b1, b2 = float(_huber(y1, Xv)[j]), float(_huber(y2, Xv)[j])
            lo, hi, fl = _boot(lambda i: _huber(y1[i], Xv[i])[j] - _huber(y2[i], Xv[i])[j], len(y1), rng)
            out["age_slope"] = dict(n=len(q), slope_fit_1_40=b1, slope_fit_2_40=b2,
                                    std_slope_fit_1_40=b1 * q.age.std(ddof=1) / q[e1].std(ddof=1),
                                    std_slope_fit_2_40=b2 * q.age.std(ddof=1) / q[e2].std(ddof=1),
                                    difference=b1 - b2, diff_ci=[lo, hi], diff_ci_excludes_zero=bool(lo > 0 or hi < 0),
                                    boot_failures=fl)
            # attention coefficient: the pre-specified fallback model, coefficient by name
            q = d[[e1, e2, "attention", "age", "sex"]].dropna().copy()
            zf = lambda v: (v - v.mean()) / v.std(ddof=1)
            q["_a"] = zf(q.attention)
            X = _design(q, f"1 + _a + {SPLINE} + C(sex)"); j = list(X.columns).index("_a"); Xv = X.values
            y1, y2 = zf(q[e1]).values, zf(q[e2]).values
            c1, c2 = float(_huber(y1, Xv)[j]), float(_huber(y2, Xv)[j])
            lo, hi, fl = _boot(lambda i: _huber(y1[i], Xv[i])[j] - _huber(y2[i], Xv[i])[j], len(y1), rng)
            out["attention_coef"] = dict(n=len(q), units="SD of exponent per SD of attention",
                                         coef_fit_1_40=c1, coef_fit_2_40=c2, difference=c1 - c2,
                                         diff_ci=[lo, hi], diff_ci_excludes_zero=bool(lo > 0 or hi < 0),
                                         boot_failures=fl)
            R2["versions"][ver] = out
        prim = R2["versions"]["matched_primary"]
        R2["decision_rule"] = dict(
            rule="claim 'fit range changes conclusions' retained only if, under the matched-set primary, the paired "
                 "95% CI for the change in age slope OR in attention coefficient excludes zero",
            age_slope_ci_excludes_zero=prim["age_slope"]["diff_ci_excludes_zero"],
            attention_ci_excludes_zero=prim["attention_coef"]["diff_ci_excludes_zero"])
        R2["decision_rule"]["claim_retained"] = bool(R2["decision_rule"]["age_slope_ci_excludes_zero"] or
                                                     R2["decision_rule"]["attention_ci_excludes_zero"])
        json.dump(R2, open(os.path.join(RUN, "P0-2_results.json"), "w"), indent=2, default=float)

        # ------------------------------------------------------------ P0-4
        Pp = pd.read_csv(paths["participants"], usecols=["participant_id", "sfreq", "n_samples"])
        r3 = CH[CH.invalid_reasons.fillna("").str.contains("stuck")].merge(Pp, on="participant_id", how="left")
        r3 = pd.DataFrame(dict(release="R3", participant_id=r3.participant_id, channel=r3.channel,
                               run_value=r3.run_value.astype(float),
                               whole=r3.longest_identical_run_s >= r3.n_samples / r3.sfreq - 1.0 / r3.sfreq))
        sc = pd.concat([pd.read_csv(f) for f in scan_files], ignore_index=True) if scan_files else pd.DataFrame()
        if len(sc):
            sc = sc[_bool(sc.stuck)]
            sc = pd.DataFrame(dict(release=sc.release, participant_id=sc.participant_id, channel=sc.channel,
                                   run_value=sc.run_value.astype(float), whole=_bool(sc.whole_recording)))
        F4 = pd.concat([r3, sc], ignore_index=True)
        F4["peripheral"] = F4.channel.isin(PERIPH)
        vals = set(F4.run_value.tolist()); rows = []
        for v, g in F4.groupby("run_value"):
            f32 = np.float32(v)
            rows.append(dict(value_repr=repr(float(v)), value_17sig=f"{v:.17g}",
                             float32_hex=struct.pack(">f", f32).hex(), exactly_float32=bool(float(f32) == v),
                             mirror_present=bool((-v) in vals and v != 0),
                             channels=len(g), children=int(g[["release", "participant_id"]].drop_duplicates().shape[0]),
                             releases=",".join(sorted(set(g.release))),
                             by_release=";".join(f"{r}:{n}" for r, n in g.release.value_counts().sort_index().items()),
                             whole=int(g.whole.sum()), episodic=int((~g.whole).sum()),
                             scalp=int((~g.peripheral).sum()), peripheral=int(g.peripheral.sum())))
        V4 = pd.DataFrame(rows).sort_values("channels", ascending=False)
        V4.to_csv(os.path.join(RUN, "P0-4_stuck_values.csv"), index=False)
        r3_220 = F4[(F4.release == "R3") & (np.abs(np.abs(F4.run_value) - 220853.4) < 1.0)]
        log["P0-4_R3_values_near_220853"] = {repr(float(v)): int(n) for v, n in r3_220.run_value.value_counts().items()}
        log["status"] = "complete"
    except GateFailure as e:
        log["status"] = f"STOPPED: {e}"
        write_log(); print(f"\nRUN STOPPED AT A GATE. Nothing further computed.\n{e}\nLog: {RUN}/run_log.json")
        raise
    write_log()

    # ---------------------------------------------------------------- printed summary (read back from files)
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70)
    print(f"RUN FOLDER: {RUN}\n")
    lg = json.load(open(os.path.join(RUN, "run_log.json")))
    print("GATES (run_log.json):")
    print(f"  G1 reproduction: all mismatches 0; max abs diff "
          f"{max(v['max_abs_diff'] for v in lg['gates']['G1'].values()):.2e}")
    print(f"  G2 C+H+G = D: max abs diff {lg['gates']['G2']['max_abs_diff']:.2e}")
    print(f"  G3 D = G when nothing excluded: n = {lg['gates']['G3']['n_nothing_excluded']}, "
          f"max abs diff {lg['gates']['G3']['max_abs_diff']:.2e}")
    s1 = pd.read_csv(os.path.join(RUN, "P0-1_summary.csv"))
    print("\nP0-1 DECOMPOSITION (P0-1_summary.csv) -- primary order and Shapley:")
    print(s1[~s1.term.str.startswith("order")][["term", "n", "median", "ci_low", "ci_high", "mean", "sd",
                                                "pct_exactly_zero"]].round(5).to_string(index=False))
    print("\n  order sensitivity, median contribution by factor (P0-1_summary.csv):")
    od = s1[s1.term.str.startswith("order")].copy()
    od["order"] = od.term.str.extract(r"order (\S+):")[0]; od["factor"] = od.term.str.split(": ").str[1]
    print(od.pivot(index="order", columns="factor", values="median").round(5).to_string())
    print("\nP0-1 SUBGROUPS (P0-1_subgroups.csv):")
    print(pd.read_csv(os.path.join(RUN, "P0-1_subgroups.csv")).round(5).to_string(index=False))
    r2 = json.load(open(os.path.join(RUN, "P0-2_results.json")))
    print(f"\nP0-2 FIT RANGE (P0-2_results.json); SD units use SD = {r2['sd_units_reference']['sd']:.4f}; "
          f"empty matched sets: {r2['n_empty_matched_set']}")
    for ver, o in r2["versions"].items():
        dd, ag = o["delta"], o["agreement"]
        print(f"  [{ver}] shift 1-40 minus 2-40: median {dd['median']:+.4f}, mean {dd['mean']:+.4f}, SD {dd['sd']:.4f}, "
              f"IQR [{dd['q25']:+.4f}, {dd['q75']:+.4f}], 5-95% [{dd['p5']:+.4f}, {dd['p95']:+.4f}], n={dd['n']}")
        print(f"      SD units: median {o['delta_sd_units']['median']:+.3f}, SD {o['delta_sd_units']['sd']:.3f}")
        print(f"      agreement: Spearman {ag['spearman']['est']:.3f} {np.round(ag['spearman']['ci'], 3).tolist()}, "
              f"ICC(3,1) {ag['icc_3_1']['est']:.3f} {np.round(ag['icc_3_1']['ci'], 3).tolist()}, "
              f"ICC(2,1) {ag['icc_2_1']['est']:.3f} {np.round(ag['icc_2_1']['ci'], 3).tolist()}")
        for p, v in o["delta_vs_predictors"].items():
            print(f"      shift ~ {p:<12} rho {v['spearman']:+.3f} {np.round(v['ci'], 3).tolist()} n={v['n']}")
        a, t = o["age_slope"], o["attention_coef"]
        print(f"      age slope /yr: 1-40 {a['slope_fit_1_40']:+.5f}, 2-40 {a['slope_fit_2_40']:+.5f}, "
              f"diff {a['difference']:+.5f} {np.round(a['diff_ci'], 5).tolist()} n={a['n']}")
        print(f"      attention coef: 1-40 {t['coef_fit_1_40']:+.3f}, 2-40 {t['coef_fit_2_40']:+.3f}, "
              f"diff {t['difference']:+.3f} {np.round(t['diff_ci'], 3).tolist()} n={t['n']}")
    dr = r2["decision_rule"]
    print(f"  DECISION RULE: age CI excludes 0 = {dr['age_slope_ci_excludes_zero']}, attention CI excludes 0 = "
          f"{dr['attention_ci_excludes_zero']}  ->  claim retained = {dr['claim_retained']}")
    print("\nP0-4 EXACT FROZEN VALUES (P0-4_stuck_values.csv):")
    print(pd.read_csv(os.path.join(RUN, "P0-4_stuck_values.csv"), dtype={"float32_hex": str}).to_string(index=False))
    print(f"\n  Release 3 values with |value| near 220,853 (run_log.json): {lg['P0-4_R3_values_near_220853']}")


if __name__ == "__main__":
    main()
