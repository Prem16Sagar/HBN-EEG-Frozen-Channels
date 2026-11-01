# =============================================================================
# HBN Release 3 -- STREAMING RUNNER for pipeline v3.1   (Colab cell)
#
# Processes participants ONE AT A TIME: fetch 4 small files -> pilot_v31.process()
# -> write per-participant results to Google Drive -> delete the raw files.
# Never needs more than ~110 MB of local disk.
#
# CRASH-SAFE + RESUMABLE: each participant gets its own folder; its participant.csv is
# written LAST (via atomic rename) and marks completion. After a Colab restart, re-run
# this cell: completed participants are skipped; nothing is duplicated.
#
# Before running:
#   1. Mount Drive:  from google.colab import drive; drive.mount('/content/drive')
#   2. Put pilot_v31.py in CODE_DIR (below).
#   3. First run with MAX_NEW = 3 to check, then set MAX_NEW = None for the full cohort.
# =============================================================================
import os, sys, time, json, shutil, tempfile, traceback, gc
import numpy as np, pandas as pd, requests

CODE_DIR = "/content/drive/MyDrive/hbn"           # folder containing pilot_v31.py
OUT      = "/content/drive/MyDrive/hbn/R3_v31"    # persistent results
CACHE    = None      # e.g. "/content/drive/MyDrive/hbn/raw_R3" to KEEP raw files (~100 MB each)
S3       = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG/cmi_bids_R3"
LOCAL_PARTICIPANTS_TSV = None   # fallback path if participants.tsv cannot be fetched from S3
INCLUDE_STATUS = ("available", "caution")   # curator RestingState status; 'caution' kept + recorded
MAX_NEW  = None       # participants to process THIS run (None = all remaining; e.g. 60 for chunks)
RETRY_FAILED = False  # True = retry participants that failed previously
REQUIRED = ["_eeg.set", "_events.tsv", "_channels.tsv"]   # what the pipeline actually reads
OPTIONAL = ["_eeg.json"]                                   # provenance only; absence is recorded
META_COLS = ["release_number", "sex", "age", "ehq_total", "p_factor", "attention",
             "internalizing", "externalizing", "RestingState"]


def atomic_write_csv(df, path):
    tmp = path + ".tmp"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def fetch_file(url, dest, retries=3):
    """HTTPS download with retries. Size is checked only for non-encoded responses
    (requests decompresses gzip, so encoded sizes legitimately differ)."""
    last = None
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                if r.status_code == 404:
                    raise FileNotFoundError(f"404 not found: {url}")
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                encoded = r.headers.get("content-encoding", "").lower() not in ("", "identity")
                done = 0
                with open(dest + ".part", "wb") as f:
                    for b in r.iter_content(chunk_size=1 << 20):
                        f.write(b); done += len(b)
            if total and not encoded and done != total:
                raise IOError(f"incomplete: {done} of {total} bytes")
            os.replace(dest + ".part", dest)
            return
        except FileNotFoundError:
            raise                                   # missing file: do not retry
        except Exception as e:
            last = e
            time.sleep(2 * attempt)
    raise IOError(f"failed after {retries} attempts: {last}")


def fetch_participant(pid, dest_dir):
    """Fetch one participant's RestingState files (from CACHE if present).
    Required files raise on failure; optional ones are skipped and reported."""
    missing_optional = []
    for sfx in REQUIRED + OPTIONAL:
        fname = f"sub-{pid}_task-RestingState{sfx}"
        dest = os.path.join(dest_dir, fname)
        cached = os.path.join(CACHE, f"sub-{pid}", fname) if CACHE else None
        if cached and os.path.exists(cached):
            shutil.copy2(cached, dest)
            continue
        try:
            fetch_file(f"{S3}/sub-{pid}/eeg/{fname}", dest)
        except FileNotFoundError:
            if sfx in OPTIONAL:
                missing_optional.append(sfx)
                continue
            raise
        if cached:
            os.makedirs(os.path.dirname(cached), exist_ok=True)
            shutil.copy2(dest, cached)
    return missing_optional


def load_participants():
    local = os.path.join(OUT, "participants.tsv")
    if not os.path.exists(local):
        if LOCAL_PARTICIPANTS_TSV:
            shutil.copy2(LOCAL_PARTICIPANTS_TSV, local)
        else:
            fetch_file(f"{S3}/participants.tsv", local)
    p = pd.read_csv(local, sep="\t")
    p["pid"] = p["participant_id"].str.replace("^sub-", "", regex=True)
    return p


def combine():
    """Merge all per-participant files into cohort-level CSVs (safe to call any time)."""
    base = os.path.join(OUT, "per_participant")
    P, C, B, F = [], [], [], []
    for pid in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        d = os.path.join(base, pid)
        if os.path.exists(os.path.join(d, "participant.csv")):
            P.append(pd.read_csv(os.path.join(d, "participant.csv")))
            C.append(pd.read_csv(os.path.join(d, "channels.csv")))
            B.append(pd.read_csv(os.path.join(d, "blocks.csv")))
        elif os.path.exists(os.path.join(d, "failed.json")):
            with open(os.path.join(d, "failed.json")) as fh:
                F.append(json.load(fh))
    if P:
        atomic_write_csv(pd.concat(P, ignore_index=True), os.path.join(OUT, "R3_v31_participants.csv"))
        atomic_write_csv(pd.concat(C, ignore_index=True), os.path.join(OUT, "R3_v31_channels.csv"))
        atomic_write_csv(pd.concat(B, ignore_index=True), os.path.join(OUT, "R3_v31_blocks.csv"))
    atomic_write_csv(pd.DataFrame(F, columns=["participant_id", "stage", "reason"]),
                     os.path.join(OUT, "R3_v31_failures.csv"))
    return len(P), len(F)


def main():
    if OUT.startswith("/content/drive") and not os.path.ismount("/content/drive"):
        raise RuntimeError("Google Drive is not mounted -- run drive.mount('/content/drive') first, "
                           "or results will be lost on restart.")
    sys.path.insert(0, CODE_DIR)
    if not os.path.exists(os.path.join(CODE_DIR, "pilot_v31.py")):
        raise FileNotFoundError(f"pilot_v31.py not found in {CODE_DIR}")
    import pilot_v31 as v
    base = os.path.join(OUT, "per_participant"); os.makedirs(base, exist_ok=True)

    parts = load_participants()
    todo = parts[parts["RestingState"].isin(INCLUDE_STATUS)]
    print(f"participants.tsv: {len(parts)} rows | RestingState in {INCLUDE_STATUS}: {len(todo)}")

    def state(pid):
        d = os.path.join(base, pid)
        if os.path.exists(os.path.join(d, "participant.csv")): return "done"
        if os.path.exists(os.path.join(d, "failed.json")): return "failed"
        return "pending"
    states = {pid: state(pid) for pid in todo["pid"]}
    queue = [pid for pid in todo["pid"] if states[pid] == "pending" or (RETRY_FAILED and states[pid] == "failed")]
    print(f"already done: {sum(s == 'done' for s in states.values())} | failed: "
          f"{sum(s == 'failed' for s in states.values())} | to process now: "
          f"{len(queue) if MAX_NEW is None else min(MAX_NEW, len(queue))} of {len(queue)} remaining\n")
    if MAX_NEW is not None:
        queue = queue[:MAX_NEW]

    times = []
    for k, pid in enumerate(queue, 1):
        d = os.path.join(base, pid); os.makedirs(d, exist_ok=True)
        if os.path.exists(os.path.join(d, "failed.json")):
            os.remove(os.path.join(d, "failed.json"))          # retrying
        tmp = tempfile.mkdtemp(prefix=f"{pid}_", dir="/content" if os.path.isdir("/content") else None)
        t0, stage = time.time(), "fetch"
        try:
            missing_optional = fetch_participant(pid, tmp)
            stage = "process"
            row, CH, BL = v.process(pid, tmp)
            meta = todo[todo["pid"] == pid].iloc[0]
            for c in META_COLS:
                row[c] = meta.get(c, np.nan)
            row["curator_resting_status"] = meta.get("RestingState", np.nan)
            row["missing_optional_files"] = ";".join(missing_optional)
            row.update(excluded=False, exclusion_reason="", pipeline_version=v.PARAMS["version"])
            stage = "write"
            atomic_write_csv(CH, os.path.join(d, "channels.csv"))
            atomic_write_csv(BL, os.path.join(d, "blocks.csv"))
            atomic_write_csv(pd.DataFrame([row]), os.path.join(d, "participant.csv"))   # LAST = done
            msg = (f"OK   exp={row['scalp__ap2_exponent']:.3f} invalid={row['n_invalid']} "
                   f"ref={row['n_reference_channels']}")
        except Exception as e:
            with open(os.path.join(d, "failed.json"), "w") as fh:
                json.dump({"participant_id": pid, "stage": stage, "reason": f"{type(e).__name__}: {e}"}, fh)
            msg = f"FAIL at {stage}: {type(e).__name__}: {e}"
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
            gc.collect()
        times.append(time.time() - t0)
        eta = np.mean(times) * (len(queue) - k) / 60
        print(f"[{k}/{len(queue)}] {pid}  {msg}  ({times[-1]:.0f} s, ETA {eta:.0f} min)")

    n_ok, n_fail = combine()
    with open(os.path.join(OUT, "R3_v31_provenance.json"), "w") as fh:
        json.dump(dict(params=v.PARAMS, groups=v.GROUPS, frontal=sorted(v.FRONTAL),
                       include_status=INCLUDE_STATUS, s3=S3, fooof=v.fooof.__version__,
                       numpy=np.__version__, pandas=pd.__version__), fh, indent=2, default=str)
    print(f"\nCohort so far: {n_ok} completed, {n_fail} failed. Combined CSVs written to {OUT}/")


if __name__ == "__main__":
    main()
