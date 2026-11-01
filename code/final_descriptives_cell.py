# =============================================================================
# FINAL DESCRIPTIVES (items 1 and 2)  -- read-only, ~1 minute
#  1. Stuck-value TYPES in the Release 3 cohort, classified exactly as in the
#     cross-release scan, so the paper's failure-type table covers all releases.
#  2. 95% bootstrap CIs for EACH attention coefficient (standard and uncensored
#     exponent). Descriptive addition: the pre-specified DIFFERENCE CI in
#     results.json remains the reported one (it came from one continuing random
#     stream; a fresh seed can move it in the third decimal).
# =============================================================================
import os, json, warnings
import numpy as np, pandas as pd, patsy, statsmodels.api as sm

R3 = "/content/drive/MyDrive/hbn/R3_v31"
OUT = os.path.join(R3, "analysis")
KNOWN, SAT = 193893.78125, 220853.4
N_BOOT, SEED = 2000, 20261020
SPLINE = "cr(age, df=3, constraints='center')"

# ---------------------------------------------------------------- item 1
# (never name a variable C: patsy formulas below use the C() function)
chans = pd.read_csv(os.path.join(R3, "R3_v31_channels.csv"),
                usecols=["participant_id", "channel", "location", "invalid_reasons",
                         "run_value", "longest_identical_run_s"])
P = pd.read_csv(os.path.join(R3, "R3_v31_participants.csv"), usecols=["participant_id", "sfreq", "n_samples"])
chans = chans.merge(P, on="participant_id", how="left")
S = chans[chans.invalid_reasons.fillna("").str.contains("stuck")].copy()
S["whole_recording"] = S.longest_identical_run_s >= S.n_samples / S.sfreq - 1.0 / S.sfreq
S["type"] = np.where(np.isclose(S.run_value, KNOWN, atol=0.05), "193,893.78 code",
            np.where(np.isclose(S.run_value.abs(), SAT, atol=1.0), "±220,853 (saturation-like)", "other constant"))
T = S.groupby("type").agg(channels=("channel", "size"), children=("participant_id", "nunique"),
                          whole_recording=("whole_recording", "sum"),
                          scalp=("location", lambda l: int((l == "scalp").sum())),
                          peripheral=("location", lambda l: int((l == "peripheral").sum())),
                          commonest=("channel", lambda c: ", ".join(f"{k} ({v})" for k, v in c.value_counts().head(4).items())))
T["episodic"] = T.channels - T.whole_recording
n_children = chans.participant_id.nunique()
print(f"1. RELEASE 3 STUCK-VALUE TYPES  (cohort n = {n_children})")
print(T[["channels", "children", "whole_recording", "episodic", "scalp", "peripheral", "commonest"]].to_string())
others = S[S.type == "other constant"].run_value.round(2).value_counts().head(8)
if len(others):
    print("   other stuck values: " + ", ".join(f"{v:g} ({k})" for v, k in others.items()))
print(f"   children with ANY stuck channel: {S.participant_id.nunique()} of {n_children}")
T.reset_index().to_csv(os.path.join(OUT, "R3_stuck_value_types.csv"), index=False)

# ---------------------------------------------------------------- item 2
A = pd.read_csv(os.path.join(OUT, "analysis_sample.csv"))
d = A[["scalp__ap2_exponent", "scalp_allvalid__ap2_exponent_nofilter_mean", "attention", "age", "sex"]].dropna().copy()
z = lambda x: (x - x.mean()) / x.std(ddof=1)
d["_a"] = z(d["attention"])
X = patsy.dmatrix("1 + _a + " + SPLINE + " + C(sex)", d, return_type="dataframe",
                  eval_env=patsy.EvalEnvironment([{}]))   # immune to names left in the notebook
j = list(X.columns).index("_a"); Xv = X.values
ys, yu = z(d["scalp__ap2_exponent"]).values, z(d["scalp_allvalid__ap2_exponent_nofilter_mean"]).values

def beta(y, Xm):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sm.RLM(y, Xm, M=sm.robust.norms.HuberT()).fit().params[j]

bs, bu = beta(ys, Xv), beta(yu, Xv)
rng = np.random.default_rng(SEED); n = len(ys); B = []
for _ in range(N_BOOT):
    i = rng.integers(0, n, n)
    try:
        B.append((beta(ys[i], Xv[i]), beta(yu[i], Xv[i])))
    except Exception:
        pass
B = np.array(B)
ci = lambda v: np.percentile(v, [2.5, 97.5])
with open(os.path.join(OUT, "results.json")) as fh:
    pre = json.load(fh)["fallback"]["secondary_attention"]
print(f"\n2. ATTENTION COEFFICIENTS (n = {n}; {len(B)} bootstrap resamples)")
print(f"   standard exponent    beta = {bs:+.3f}  95% CI [{ci(B[:, 0])[0]:+.3f}, {ci(B[:, 0])[1]:+.3f}]")
print(f"   uncensored exponent  beta = {bu:+.3f}  95% CI [{ci(B[:, 1])[0]:+.3f}, {ci(B[:, 1])[1]:+.3f}]")
print(f"   difference (this run) {bs - bu:+.3f}  95% CI [{ci(B[:, 0] - B[:, 1])[0]:+.3f}, {ci(B[:, 0] - B[:, 1])[1]:+.3f}]")
print(f"   difference (PRE-SPECIFIED, reported) {pre['difference']:+.3f}  95% CI [{pre['ci_low']:+.3f}, {pre['ci_high']:+.3f}]")
with open(os.path.join(OUT, "attention_individual_cis.json"), "w") as fh:
    json.dump(dict(n=n, n_boot=len(B), seed=SEED, beta_standard=bs, ci_standard=ci(B[:, 0]).tolist(),
                   beta_uncensored=bu, ci_uncensored=ci(B[:, 1]).tolist(),
                   note="descriptive addition; the pre-specified difference CI in results.json is the reported one"),
              fh, indent=2)
print(f"\nSaved: {OUT}/R3_stuck_value_types.csv and attention_individual_cis.json")
