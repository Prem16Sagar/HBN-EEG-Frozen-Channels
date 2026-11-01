# =============================================================================
# AMENDMENT A2, STEP 2 -- RELEASES 10 & 11 FROM OPENNEURO      (Colab cell)
# Implements revision_B_A2_step2_spec_v2.md (approved).
#
# * The census module wp1_census.py is IMPORTED UNCHANGED from CODE_DIR; only its
#   download addresses are redirected (FCP-INDI cmi_bids_Rn  ->  OpenNeuro dsXXXXXX).
# * Gates, before any Release 10/11 scanning (failure stops the run):
#     G-C3  sidecar files (10 children/release) + one .set per release identical on
#           FCP-INDI and OpenNeuro, Releases 1-9
#     G-C4  3 Release 3 children re-scanned from OpenNeuro == census_v1 results exactly
# * Population: every sub-* folder on OpenNeuro for ds005515/ds005516; a folder without a
#   resting-state .set is recorded as "no resting-state file" (not a failure).
# * Combined summary R1-R11 in a NEW folder; releases with different selection rules are
#   reported separately, never pooled. Gate G-C5: R1-R9 counts must reproduce census_v1.
# Run after: drive.mount('/content/drive'). wp1_census.py must be in CODE_DIR.
# =============================================================================
import os, sys, json, time, glob, hashlib, xml.etree.ElementTree as ET
import numpy as np, pandas as pd, requests
from io import StringIO
from scipy import stats

DRIVE    = "/content/drive/MyDrive/hbn"
CODE_DIR = DRIVE
CENSUS   = f"{DRIVE}/revision_B/census/census_v1"
RUN      = f"{DRIVE}/revision_B/census/census_v1_r10r11"
FCP      = "https://fcp-indi.s3.amazonaws.com"
FCP_PFX  = "data/Projects/HBN/BIDS_EEG"
ON       = "https://s3.amazonaws.com/openneuro.org"
DS = {"R1": "ds005505", "R2": "ds005506", "R3": "ds005507", "R4": "ds005508", "R5": "ds005509",
      "R6": "ds005510", "R7": "ds005511", "R8": "ds005512", "R9": "ds005514",
      "R10": "ds005515", "R11": "ds005516"}
NEW = ("R10", "R11")
G_C3_N, G_C4_N, MAX_NEW = 10, 3, None
SENTINEL = "SENTINEL://fcp"
NSX = "{http://s3.amazonaws.com/doc/2006-03-01/}"


class GateFailure(RuntimeError):
    pass


def http_get(url, **kw):
    return requests.get(url, timeout=120, **kw)


def sha256_url(url):
    r = http_get(url, stream=True)
    if r.status_code != 200:
        return dict(status=r.status_code)
    h, n = hashlib.sha256(), 0
    for b in r.iter_content(1 << 20):
        h.update(b); n += len(b)
    return dict(status=200, bytes=n, sha256=h.hexdigest())


def list_subjects(ds):
    subs, token = [], None
    while True:
        url = f"{ON}/?list-type=2&prefix={ds}/&delimiter=/" + (f"&continuation-token={requests.utils.quote(token)}" if token else "")
        r = http_get(url)
        if r.status_code != 200:
            raise RuntimeError(f"listing {ds} failed: HTTP {r.status_code}")
        root = ET.fromstring(r.content)
        for p in root.findall(f"{NSX}CommonPrefixes"):
            name = p.find(f"{NSX}Prefix").text.rstrip("/").split("/")[-1]
            if name.startswith("sub-"): subs.append(name[4:])
        nxt = root.find(f"{NSX}NextContinuationToken")
        if nxt is None: return sorted(subs)
        token = nxt.text


def load_census_module():
    sys.path.insert(0, CODE_DIR)
    import wp1_census as c
    original_fetch = c.fetch
    def redirected(url, dest, retries=3):
        for rel, ds in DS.items():
            pre = f"{SENTINEL}/cmi_bids_{rel}/"
            if url.startswith(pre):
                url = f"{ON}/{ds}/" + url[len(pre):]; break
        return original_fetch(url, dest, retries)
    c.S3_ROOT, c.fetch = SENTINEL, redirected
    return c


def atomic_json(obj, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, default=float)
    os.replace(path + ".tmp", path)


def census_children():
    d = os.path.join(CENSUS, "per_child"); out = {}
    for f in sorted(os.listdir(d)):
        if f.endswith(".json") and not f.endswith(".failed.json"):
            rel, pid = f[:-5].split("_", 1); out.setdefault(rel, []).append(pid)
    return out


def gate_c3(cc):
    rows = []
    for rel in [f"R{i}" for i in range(1, 10)]:
        pids = sorted(cc.get(rel, []))[:G_C3_N]
        for k, pid in enumerate(pids):
            sfx = ["_channels.tsv", "_events.tsv"] + (["_eeg.set"] if k == 0 else [])
            for s in sfx:
                rp = f"sub-{pid}/eeg/sub-{pid}_task-RestingState{s}"
                a = sha256_url(f"{FCP}/{FCP_PFX}/cmi_bids_{rel}/{rp}"); b = sha256_url(f"{ON}/{DS[rel]}/{rp}")
                rows.append(dict(release=rel, participant_id=pid, file=s, fcp_status=a.get("status"),
                                 on_status=b.get("status"), identical=a.get("sha256") is not None and a.get("sha256") == b.get("sha256")))
    T = pd.DataFrame(rows)
    return T, bool(len(T)) and bool(T.identical.all())


def gate_c4(c, cc):
    res = []
    for pid in sorted(cc.get("R3", []))[:G_C4_N]:
        with open(os.path.join(CENSUS, "per_child", f"R3_{pid}.json")) as fh: old = json.load(fh)
        new = json.loads(json.dumps(c.scan_child("R3", pid), default=float))      # same JSON round-trip as storage
        res.append(dict(participant_id=pid, identical=(old == new)))
    return res, bool(res) and all(r["identical"] for r in res)


def load_records():
    recs = []
    for folder in (os.path.join(CENSUS, "per_child"), os.path.join(RUN, "per_child")):
        if not os.path.isdir(folder): continue
        for f in sorted(os.listdir(folder)):
            if f.endswith(".json") and not (f.endswith(".failed.json") or f.endswith(".nofile.json")):
                with open(os.path.join(folder, f)) as fh: recs.append(json.load(fh))
    return recs


def flags(j, KNOWN):
    fz = j["frozen"]; sc = [c for c in fz if c["location"] == "scalp"]
    return dict(release=j["release"], participant_id=j["participant_id"],
                code=any(float(v) == KNOWN for c in fz for v in c["values"]), any_frozen=len(fz) > 0,
                episodic=any(not c["whole_recording"] for c in fz),
                scalp_any=len(sc) > 0, scalp_episodic=any(not c["whole_recording"] for c in sc),
                scalp_whole_only=len(sc) > 0 and all(c["whole_recording"] for c in sc))


def main():
    if RUN.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Mount Google Drive first: drive.mount('/content/drive')")
    if not os.path.exists(os.path.join(CODE_DIR, "wp1_census.py")):
        raise FileNotFoundError(f"put wp1_census.py in {CODE_DIR} (it is imported unchanged)")
    c = load_census_module()
    os.makedirs(os.path.join(RUN, "per_child"), exist_ok=True)
    cc = census_children()

    # ---- provenance of the OpenNeuro copies (version is not pinned by S3)
    prov = {}
    for rel in NEW:
        for fn in ("dataset_description.json", "participants.tsv", "CHANGES"):
            prov[f"{rel}/{fn}"] = sha256_url(f"{ON}/{DS[rel]}/{fn}")
    # ---- gates G-C3, G-C4
    T3, ok3 = gate_c3(cc)
    g4, ok4 = gate_c4(c, cc)
    S = os.path.join(RUN, f"summary_{time.strftime('%Y%m%d_%H%M%S')}"); os.makedirs(S, exist_ok=False)
    T3.to_csv(os.path.join(S, "gate_G_C3.csv"), index=False)
    gates = {"G-C3": dict(files=len(T3), identical=int(T3.identical.sum()) if len(T3) else 0, passed=ok3),
             "G-C4": dict(children=len(g4), detail=g4, passed=ok4)}
    atomic_json(dict(gates=gates, provenance=prov), os.path.join(S, "gates_provenance.json"))
    print(f"GATE G-C3: {gates['G-C3']}\nGATE G-C4: {[(g['participant_id'], g['identical']) for g in g4]} passed={ok4}")
    if not (ok3 and ok4):
        raise GateFailure(f"gate failed -- nothing scanned; see {S}")

    # ---- scan Releases 10 and 11 (every subject folder)
    d = os.path.join(RUN, "per_child"); todo = []
    for rel in NEW:
        for pid in list_subjects(DS[rel]):
            done = any(os.path.exists(os.path.join(d, f"{rel}_{pid}{e}")) for e in (".json", ".nofile.json", ".failed.json"))
            if not done: todo.append((rel, pid))
    print(f"{len(todo)} subject folders to process")
    if MAX_NEW is not None: todo = todo[:MAX_NEW]
    t_all = time.time()
    for k, (rel, pid) in enumerate(todo, 1):
        try:
            atomic_json(c.scan_child(rel, pid), os.path.join(d, f"{rel}_{pid}.json")); msg = "OK"
        except FileNotFoundError:
            atomic_json(dict(release=rel, participant_id=pid, reason="no resting-state file"), os.path.join(d, f"{rel}_{pid}.nofile.json")); msg = "no resting file"
        except Exception as e:
            atomic_json(dict(release=rel, participant_id=pid, reason=f"{type(e).__name__}: {e}"), os.path.join(d, f"{rel}_{pid}.failed.json")); msg = f"FAIL {e}"
        if k % 25 == 0 or k == len(todo) or msg.startswith("FAIL"):
            print(f"  [{k}/{len(todo)}] {rel} {pid} {msg} (ETA {(time.time() - t_all) / k * (len(todo) - k) / 60:.0f} min)")

    # ---- combined summary
    F = pd.DataFrame([flags(j, c.KNOWN) for j in load_records()])
    nofile = [f for f in os.listdir(d) if f.endswith(".nofile.json")]
    failed = [f for f in os.listdir(d) if f.endswith(".failed.json")]
    rule = {r: "census rule (participants.tsv available/caution)" for r in [f"R{i}" for i in range(1, 10)]}
    rule.update({r: "folder rule (every subject folder with a resting-state file)" for r in NEW})
    rows = []
    def add(label, g, rl):
        n = len(g); row = dict(release=label, selection_rule=rl, scanned=n)
        for col in ("code", "any_frozen", "episodic", "scalp_any", "scalp_episodic", "scalp_whole_only"):
            k_ = int(g[col].sum()); lo, hi = c.wilson(k_, n)
            row.update({f"n_{col}": k_, f"pct_{col}": 100 * k_ / n if n else np.nan,
                        f"ci_{col}": f"[{max(0.0, 100 * lo):.1f}, {min(100.0, 100 * hi):.1f}]" if n else ""})
        rows.append(row)
    for rel in [f"R{i}" for i in range(1, 12)]:
        add(rel, F[F.release == rel], rule[rel])
    # Release 10 under the census rule as well (cross-check)
    p10 = pd.read_csv(StringIO(http_get(f"{ON}/{DS['R10']}/participants.tsv").text), sep="\t")
    av = set(p10[p10.RestingState.isin(["available", "caution"])].participant_id.str.replace("^sub-", "", regex=True))
    add("R10 (census rule)", F[(F.release == "R10") & F.participant_id.isin(av)], rule["R1"])
    add("R1-R9 total", F[F.release.isin([f"R{i}" for i in range(1, 10)])], "census rule")
    T = pd.DataFrame(rows); T.to_csv(os.path.join(S, "combined_summary.csv"), index=False)
    scanned10 = set(F[F.release == "R10"].participant_id)
    nofile10 = {f[4:-12] for f in nofile if f.startswith("R10_")}
    xc = dict(r10_folders_with_file=len(scanned10), r10_with_file_and_available_or_caution=len(scanned10 & av),
              r10_with_file_not_available_or_caution=len(scanned10 - av),
              r10_available_or_caution_without_file=len(av - scanned10),
              r10_folders_without_file=len(nofile10),
              r11_folders_with_file=int((F.release == "R11").sum()),
              r11_folders_without_file=sum(f.startswith("R11_") for f in nofile), failed=len(failed))
    # G-C5: R1-R9 counts must reproduce census_v1's summary
    cs = sorted(glob.glob(f"{CENSUS}/summary_*/census_summary.csv"))[-1]
    old = pd.read_csv(cs).set_index("release")
    mism = [rel for rel in [f"R{i}" for i in range(1, 10)]
            if (int(old.loc[rel, "scanned"]), int(old.loc[rel, "n_known_code"]), int(old.loc[rel, "n_any_frozen"]), int(old.loc[rel, "n_episodic"]))
            != tuple(int(T.set_index("release").loc[rel, k]) for k in ("scanned", "n_code", "n_any_frozen", "n_episodic"))]
    gates["G-C5"] = dict(compared_with=cs, mismatched_releases=mism, passed=not mism)
    chi = {}
    for label, rels in (("R1-R9", [f"R{i}" for i in range(1, 10)]), ("R1-R10 (R10 census rule)", [f"R{i}" for i in range(1, 10)] + ["R10 (census rule)"])):
        sub = T.set_index("release").loc[rels]
        sub = sub[sub["scanned"] > 0]                                    # releases with nothing scanned are excluded
        for col in ("code", "any_frozen", "scalp_episodic"):
            tab = np.array([sub[f"n_{col}"], sub["scanned"] - sub[f"n_{col}"]])
            if tab.shape[1] < 2 or (tab.sum(axis=1) == 0).any():
                chi[f"{label} | {col}"] = dict(chi2=None, dof=None, p=None, note="not computable"); continue
            x2, p, dof, _ = stats.chi2_contingency(tab)
            chi[f"{label} | {col}"] = dict(chi2=float(x2), dof=int(dof), p=float(p), releases=list(sub.index))
    atomic_json(dict(gates=gates, provenance=prov, r10_rule_crosscheck=xc, chi2=chi), os.path.join(S, "gates_provenance.json"))
    if not gates["G-C5"]["passed"]:
        raise GateFailure(f"G-C5: combined summary does not reproduce census_v1 for {mism}")

    pd.set_option("display.width", 260)
    print(f"\nSUMMARY FOLDER: {S}")
    with open(os.path.join(S, "gates_provenance.json")) as fh: G = json.load(fh)
    print("GATES (gates_provenance.json):"); [print(f"  {k}: passed={v['passed']}") for k, v in G["gates"].items()]
    print(f"\nR10/R11 FOLDERS (gates_provenance.json): {G['r10_rule_crosscheck']}")
    T = pd.read_csv(os.path.join(S, "combined_summary.csv"))
    cols = ["release", "scanned"] + [f"{p}_{k}" for k in ("code", "any_frozen", "episodic", "scalp_episodic", "scalp_whole_only") for p in ("pct", "ci")]
    print("\nCOMBINED SUMMARY (combined_summary.csv):"); print(T[cols].round(1).to_string(index=False))
    print("\nselection rules (combined_summary.csv):"); print(T[["release", "selection_rule"]].to_string(index=False))
    print("\nCHI-SQUARE (gates_provenance.json):"); [print(f"  {k}: " + (f"chi2({v['dof']}) = {v['chi2']:.1f}, p = {v['p']:.2g}" if v["chi2"] is not None else v["note"])) for k, v in G["chi2"].items()]


if __name__ == "__main__":
    main()
