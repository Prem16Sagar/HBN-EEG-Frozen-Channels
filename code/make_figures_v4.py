# =============================================================================
# FIGURES v4 for the frozen-channel paper (Option B)          (Colab cell)
# Reads only saved result files; writes PNG + PDF to a NEW folder, figures_v4.
#   Fig 1  (a) example recording's stored values; (b) frozen-channel prevalence by
#          release, Releases 1-10 + census total, Wilson 95% CIs   [combined_summary.csv]
#   Fig 2  per-child inflation and exponent error, frozen channels added to the
#          reference, by freeze type                                  [harm_per_child.csv]
#   Fig 3  PREP detection by longest frozen episode, cross-release sample, with and
#          without the one failed recording                           [detector_xrel_per_channel.csv]
# Each figure is built independently; a failure is reported and the rest still run.
# Every number drawn is also written to figure_numbers.json for the captions.
# Run after: drive.mount('/content/drive')
# =============================================================================
import os, sys, json, glob, warnings
import numpy as np, pandas as pd, scipy.io, requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore", category=RuntimeWarning)

DRIVE   = "/content/drive/MyDrive/hbn"
FIG_DIR = f"{DRIVE}/figures_v4"
EXAMPLE_PID = "NDARAA948VFH"
RAW_SEARCH  = ["/content/hbn_one", f"{DRIVE}/hbn_one", f"{DRIVE}/raw_R3/sub-{EXAMPLE_PID}"]
S3_R3   = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG/cmi_bids_R3"

# ---- style: Okabe-Ito colour-blind-safe palette, Nature-portfolio widths -------------
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7", grey="#8C8C8C")
W2 = 7.0                          # double-column width (in), <= 180 mm
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5,
                     "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300})


# ----------------------------------------------------------------------- helpers
def panel(ax, letter):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.png / .pdf")


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def longest_run(x):
    same = np.diff(x) == 0
    if not same.any():
        return 1, float(x[0]), 0
    d = np.diff(np.concatenate([[0], same.astype(np.int8), [0]]))
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    k = int(np.argmax(en - st))
    return int(en[k] - st[k] + 1), float(x[st[k]]), int(st[k])


def get_example_raw():
    """Load the example recording from a local copy, or download it once."""
    base = f"sub-{EXAMPLE_PID}_task-RestingState"
    for d in RAW_SEARCH:
        if os.path.exists(os.path.join(d, base + "_eeg.set")):
            folder = d; break
    else:
        folder = "/content/fig_example"; os.makedirs(folder, exist_ok=True)
        for sfx in ("_eeg.set", "_channels.tsv", "_events.tsv"):
            dest = os.path.join(folder, base + sfx)
            if not os.path.exists(dest):
                with requests.get(f"{S3_R3}/sub-{EXAMPLE_PID}/eeg/{base}{sfx}", stream=True, timeout=120) as r:
                    r.raise_for_status()
                    with open(dest, "wb") as f:
                        for b in r.iter_content(1 << 20): f.write(b)
    m = scipy.io.loadmat(os.path.join(folder, base + "_eeg.set"), squeeze_me=True, struct_as_record=False)
    X = np.asarray(m["data"], dtype=np.float32); sf = float(m["srate"])
    names = pd.read_csv(os.path.join(folder, base + "_channels.tsv"), sep="\t")["name"].tolist()
    ev = pd.read_csv(os.path.join(folder, base + "_events.tsv"), sep="\t")
    ev = ev[ev["value"].isin(["instructed_toCloseEyes", "instructed_toOpenEyes"])].sort_values("onset").reset_index(drop=True)
    blocks = []
    for i, r in ev.iterrows():
        if r["value"] == "instructed_toCloseEyes":
            nx = ev[(ev.index > i) & (ev["value"] == "instructed_toOpenEyes")]
            end = nx["onset"].iloc[0] if len(nx) else X.shape[1] / sf
            s, e = int((r["onset"] + 2) * sf), int((end - 2) * sf)
            if e > s: blocks.append((s, e))
    return X, sf, names, blocks




# ----------------------------------------------------------------------- inputs
def latest(pattern, must_contain=None):
    c = sorted(glob.glob(pattern))
    if must_contain:
        c = [d for d in c if os.path.exists(os.path.join(d, must_contain))]
    if not c:
        raise FileNotFoundError(f"nothing matches {pattern}" + (f" containing {must_contain}" if must_contain else ""))
    return c[-1]


def census_table():
    d = latest(f"{DRIVE}/revision_B/census/census_v1_r10r11/summary_*", "combined_summary.csv")
    T = pd.read_csv(os.path.join(d, "combined_summary.csv")).set_index("release")
    rows = [f"R{i}" for i in range(1, 10)] + ["R10 (census rule)"]
    missing = [r for r in rows if r not in T.index]
    if missing:
        raise KeyError(f"combined_summary.csv lacks rows {missing}")
    T = T.loc[rows].copy(); T.index = [r.split(" ")[0] for r in rows]
    return T, os.path.join(d, "combined_summary.csv")


# ----------------------------------------------------------------------- Figure 1
def figure1():
    from matplotlib.lines import Line2D
    T, src = census_table()
    try:
        X, sf, names, blocks = get_example_raw()
        have_raw = True
    except Exception as e:                       # raw recording unavailable: prevalence panel only
        print(f"  example recording unavailable ({type(e).__name__}: {e}); drawing panel b only")
        have_raw = False
    if have_raw:
        fig = plt.figure(figsize=(W2, 2.9))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.55], wspace=0.42)
    else:
        fig = plt.figure(figsize=(W2 * 0.62, 2.9))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.0001, 1], wspace=0)

    if have_raw:
        I = pd.DataFrame([(c,) + longest_run(X[i]) for i, c in enumerate(names)], columns=["ch", "L", "val", "start"])
        I["run_s"] = I.L / sf; dur = X.shape[1] / sf
        whole = I[(I.run_s >= dur - 1.0 / sf) & (I.ch != "Cz")].ch.tolist()
        epis = I[(I.run_s >= 1) & (I.run_s < dur - 1.0 / sf) & (I.ch != "Cz")]
        normal = next(c for c in ("E62", "E75", "E70") if c in names and c not in whole and c not in epis.ch.tolist())
        show = [(normal, "valid", OI["blue"])]
        if len(epis): show.append((epis.iloc[0].ch, "partly\nfrozen", OI["vermillion"]))
        if whole: show.append((whole[0], "frozen\nthroughout", OI["orange"]))
        if len(epis):                                   # show the whole episode where it fits (<= 20 s window)
            s0 = int(epis.iloc[0].start); e0 = s0 + int(epis.iloc[0].L)
            w0 = max(0, s0 - int(3 * sf)); w1 = min(X.shape[1], e0 + int(3 * sf), w0 + int(20 * sf))
        else:
            w0 = blocks[0][0]; w1 = min(X.shape[1], w0 + int(12 * sf))
        t = np.arange(w0, w1) / sf
        sub = gs[0, 0].subgridspec(len(show), 1, hspace=0.35)
        for k, (c, lab, col) in enumerate(show):
            ax = fig.add_subplot(sub[k])
            ax.plot(t, X[names.index(c), w0:w1].astype(float), color=col, lw=0.6)
            ax.set_ylabel(f"{c}\n{lab}", rotation=0, ha="right", va="center", fontsize=6, color=col)
            ax.tick_params(labelsize=5.5); ax.yaxis.get_offset_text().set_fontsize(5.5)
            if k < len(show) - 1: ax.set_xticklabels([])
            if k == 0: panel(ax, "a"); ax.set_title(f"Stored values, {EXAMPLE_PID}", loc="left", pad=6)
        ax.set_xlabel("time (s)")

    # (b) prevalence by release + census total, three measures side by side, Wilson CIs
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b" if have_raw else "")
    series = [("any_frozen", "any frozen channel", OI["sky"]),
              ("code", "193,893.78125 value", OI["blue"]),
              ("scalp_episodic", "partly frozen scalp channel", OI["vermillion"])]
    labels = list(T.index) + ["All"]
    N = list(T.scanned.astype(int)) + [int(T.scanned.sum())]
    xs = np.arange(len(labels)).astype(float); xs[-1] += 0.6; bw = 0.26
    for j, (col, lab, colour) in enumerate(series):
        K = list(T[f"n_{col}"].astype(int)) + [int(T[f"n_{col}"].sum())]
        p = np.array([k / n for k, n in zip(K, N)]); ci = np.array([wilson(k, n) for k, n in zip(K, N)])
        x = xs + (j - 1) * bw
        ax.bar(x, 100 * p, width=bw, color=colour, label=lab)
        ax.errorbar(x, 100 * p, yerr=[100 * (p - ci[:, 0]), 100 * (ci[:, 1] - p)], fmt="none", ecolor="k", lw=0.5, capsize=1)
    ax.axvline(xs[-2] + 0.8, color=OI["grey"], lw=0.5, ls=":")
    ax.set_xticks(xs); ax.set_xticklabels([f"{l}\n{n}" for l, n in zip(labels, N)], fontsize=5.5)
    ax.set_ylim(0, 100); ax.set_ylabel("recordings affected (%)")
    ax.set_title("Resting-state recordings with frozen channels (release / n)", loc="left")
    ax.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=3, bbox_to_anchor=(0, 1.0))
    ax.set_ylim(0, 112)
    save(fig, "Fig1_frozen_channels" if have_raw else "Fig1b_prevalence_only")
    return dict(source=src, total=int(T.scanned.sum()),
                pct={c: round(100 * T[f"n_{c}"].sum() / T.scanned.sum(), 1) for c, _, _ in series})


# ----------------------------------------------------------------------- Figure 2
GROUPS = [("no frozen scalp channel", "none"),
          ("whole-recording frozen scalp channels only", "frozen\nthroughout only"),
          (">=1 episodic frozen scalp channel", "≥1 partly\nfrozen")]


def figure2():
    d = latest(f"{DRIVE}/revision_B/harm/harm_v1/summary_*", "harm_per_child.csv")
    H = pd.read_csv(os.path.join(d, "harm_per_child.csv"))
    unknown = set(H.group) - {g for g, _ in GROUPS}
    if unknown:
        raise KeyError(f"unexpected groups in harm_per_child.csv: {unknown}")
    rng = np.random.default_rng(1)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.7), gridspec_kw=dict(wspace=0.35))
    out = {}
    for ax, col, ylab, letter in ((axs[0], "frozen_only__inflation_ratio", "amplitude inflation (×, log scale)", "a"),
                                  (axs[1], "frozen_only__exponent_error", "exponent error", "b")):
        panel(ax, letter)
        for j, (g, lab) in enumerate(GROUPS):
            v = H.loc[H.group == g, col].dropna().values
            ax.scatter(j + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=(OI["grey"], OI["orange"], OI["vermillion"])[j],
                       alpha=0.7, edgecolors="none")
            if len(v): ax.plot([j - 0.28, j + 0.28], [np.median(v)] * 2, color="k", lw=1.2)
            out[f"{col}|{g}"] = dict(n=len(v), median=float(np.median(v)) if len(v) else None)
        ax.set_xticks(range(len(GROUPS)))
        ticks = []
        for g, lab in GROUPS:
            n_all = int((H.group == g).sum()); n_ok = int(H.loc[H.group == g, col].notna().sum())
            ticks.append(f"{lab}\nn = {n_ok}" + (f"\n({n_all - n_ok} not fitted)" if n_ok < n_all else ""))
        ax.set_xticklabels(ticks, fontsize=5.5)
        ax.set_ylabel(ylab)
    axs[0].set_yscale("log"); axs[0].axhline(1, color=OI["grey"], lw=0.5, ls=":")
    axs[0].set_title("Frozen scalp channels added to the reference", loc="left")
    axs[1].axhspan(-0.05, 0.05, color=OI["grey"], alpha=0.15, lw=0)
    axs[1].axhline(0, color=OI["grey"], lw=0.5)
    axs[1].set_title("Shift in the aperiodic exponent", loc="left")
    w = H[H.group == GROUPS[1][0]].sort_values("frozen_only__inflation_ratio")
    if len(w) and w.frozen_only__inflation_ratio.iloc[-1] > 2:
        r = w.iloc[-1]
        axs[0].annotate(r.participant_id, (1, r.frozen_only__inflation_ratio), xytext=(-8, 0),
                        textcoords="offset points", fontsize=5.5, va="center", ha="right")
    save(fig, "Fig2_reference_harm")
    out["source"] = os.path.join(d, "harm_per_child.csv")
    return out


# ----------------------------------------------------------------------- Figure 3
BAND_ORDER = [("[1, 5)", "1–5 s"), ("[5, 30)", "5–30 s"), ("[30", "≥30 s"), ("whole", "whole\nrecording")]


def band_key(b):
    s = str(b).strip().lower()
    for k, (pre, _) in enumerate(BAND_ORDER):
        if s.startswith(pre.lower()):
            return k
    raise KeyError(f"unrecognised duration band {b!r}")


def figure3():
    d = latest(f"{DRIVE}/revision_B/detector_xrel/xrel_v1/summary_*", "detector_xrel_per_channel.csv")
    D = pd.read_csv(os.path.join(d, "detector_xrel_per_channel.csv"))
    D["caught"] = D["caught"].astype(str).str.lower().isin(("true", "1"))
    D["k"] = D["band"].map(band_key)
    per_rec = D.groupby("participant_id").agg(n=("caught", "size"), c=("caught", "sum"))
    worst = per_rec.assign(miss=per_rec.n - per_rec.c).sort_values("miss").index[-1]
    fig, ax = plt.subplots(figsize=(3.6, 2.6)); panel(ax, "")
    out = {}
    for j, (frame, lab, colour) in enumerate(((D, "all sampled recordings", OI["blue"]),
                                              (D[D.participant_id != worst], f"excluding {worst}", OI["green"]))):
        x = np.arange(len(BAND_ORDER)) + (j - 0.5) * 0.36
        rows = [frame[frame.k == k] for k in range(len(BAND_ORDER))]
        K = np.array([int(r.caught.sum()) for r in rows]); N = np.array([len(r) for r in rows])
        p = np.where(N > 0, K / np.maximum(N, 1), np.nan); ci = np.array([wilson(k, n) for k, n in zip(K, N)])
        ax.bar(x, 100 * p, width=0.36, color=colour, label=lab)
        ax.errorbar(x, 100 * p, yerr=[100 * (p - ci[:, 0]), 100 * (ci[:, 1] - p)], fmt="none", ecolor="k", lw=0.5, capsize=1.5)
        for xi, k, n in zip(x, K, N):
            ax.text(xi, 2, f"{k}/{n}", rotation=90, ha="center", va="bottom", fontsize=5, color="white")
        out[lab] = dict(zip([b for _, b in BAND_ORDER], zip(K.tolist(), N.tolist())))
    ax.set_xticks(range(len(BAND_ORDER))); ax.set_xticklabels([b for _, b in BAND_ORDER], fontsize=6)
    ax.set_xlabel("longest frozen episode on the channel"); ax.set_ylabel("frozen channels flagged by PREP (%)")
    ax.set_ylim(0, 115); ax.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2)
    ax.set_title("PREP detection, cross-release sample", loc="left")
    save(fig, "Fig3_detection_by_duration")
    out["worst_recording"] = str(worst); out["source"] = os.path.join(d, "detector_xrel_per_channel.csv")
    return out


def main():
    if DRIVE.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    if os.path.exists(FIG_DIR) and os.listdir(FIG_DIR):
        raise FileExistsError(f"{FIG_DIR} already has files; rename or delete it to rebuild")
    report = {}
    for name, fn in (("Figure 1", figure1), ("Figure 2", figure2), ("Figure 3", figure3)):
        print(f"{name} ...")
        try:
            report[name] = fn()
        except Exception as e:
            report[name] = f"FAILED: {type(e).__name__}: {e}"
            print(f"  {report[name]}")
    os.makedirs(FIG_DIR, exist_ok=True)
    with open(os.path.join(FIG_DIR, "figure_numbers.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(f"\nFIGURE FOLDER: {FIG_DIR}")
    print("NUMBERS DRAWN (figure_numbers.json):")
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
