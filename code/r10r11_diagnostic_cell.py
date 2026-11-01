# =============================================================================
# RELEASES 10 & 11 -- DIAGNOSTIC (amendment A2, step 1)     -- read-only, no scanning
#  1. What exists under cmi_bids_R10/ and cmi_bids_R11/ on FCP-INDI (S3 listing).
#  2. Whether participants.tsv exists for ds005515 / ds005516 on OpenNeuro.
#  3. Equivalence: SHA-256 of resting-state .set files from FCP-INDI vs OpenNeuro, for
#     three Release 3 children (ds005507) and, if present, two Release 10 children.
# Files are hashed while streaming; nothing is saved except the diagnostic report.
# Writes to revision_B/checks/r10r11_diag_<timestamp>/diagnostic.json.
# =============================================================================
import os, json, time, hashlib, xml.etree.ElementTree as ET
import pandas as pd, requests

DRIVE = "/content/drive/MyDrive/hbn"
FCP = "https://fcp-indi.s3.amazonaws.com"
FCP_PREFIX = "data/Projects/HBN/BIDS_EEG"
ON = "https://s3.amazonaws.com/openneuro.org"
DS = {"R3": "ds005507", "R10": "ds005515", "R11": "ds005516"}
OUT = f"{DRIVE}/revision_B/checks/r10r11_diag_{time.strftime('%Y%m%d_%H%M%S')}"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def http_get(url, **kw):
    return requests.get(url, timeout=120, **kw)


def list_prefix(base, prefix, max_pages=10):
    """S3 ListObjectsV2 with delimiter '/': returns (sub-folder names, root file keys)."""
    folders, files, token = [], [], None
    for _ in range(max_pages):
        url = f"{base}/?list-type=2&prefix={prefix}&delimiter=/" + (f"&continuation-token={requests.utils.quote(token)}" if token else "")
        r = http_get(url)
        if r.status_code != 200:
            return dict(status=r.status_code, folders=[], files=[])
        root = ET.fromstring(r.content)
        folders += [p.find(f"{NS}Prefix").text for p in root.findall(f"{NS}CommonPrefixes")]
        files += [c.find(f"{NS}Key").text for c in root.findall(f"{NS}Contents")]
        nxt = root.find(f"{NS}NextContinuationToken")
        if nxt is None: break
        token = nxt.text
    return dict(status=200, folders=folders, files=files)


def sha256_url(url):
    r = http_get(url, stream=True)
    if r.status_code != 200:
        return dict(status=r.status_code)
    h, n = hashlib.sha256(), 0
    for b in r.iter_content(1 << 20):
        h.update(b); n += len(b)
    return dict(status=200, bytes=n, sha256=h.hexdigest())


def set_path(pid):
    return f"sub-{pid}/eeg/sub-{pid}_task-RestingState_eeg.set"


def main():
    rep = {}
    # 1. FCP-INDI folders
    for rel in ("R10", "R11"):
        L = list_prefix(FCP, f"{FCP_PREFIX}/cmi_bids_{rel}/")
        subs = [f for f in L["folders"] if f.rstrip("/").split("/")[-1].startswith("sub-")]
        rep[f"fcp_{rel}"] = dict(listing_status=L["status"], n_subject_folders=len(subs),
                                 first_subjects=[s.rstrip("/").split("/")[-1] for s in subs[:3]],
                                 root_files=[f.split("/")[-1] for f in L["files"]],
                                 other_folders=[f.rstrip("/").split("/")[-1] for f in L["folders"] if f not in subs])
    # 2. OpenNeuro participants.tsv
    for rel in ("R10", "R11"):
        r = http_get(f"{ON}/{DS[rel]}/participants.tsv")
        info = dict(status=r.status_code)
        if r.status_code == 200:
            from io import StringIO
            t = pd.read_csv(StringIO(r.text), sep="\t")
            info.update(n_rows=len(t), has_RestingState="RestingState" in t.columns,
                        resting_available_or_caution=int(t["RestingState"].isin(["available", "caution"]).sum())
                        if "RestingState" in t.columns else None)
        rep[f"openneuro_{rel}_participants"] = info
    # 3. Equivalence: Release 3 (and Release 10 if present on FCP-INDI)
    A = pd.read_csv(f"{DRIVE}/R3_v31/analysis/analysis_sample.csv", usecols=["participant_id"])
    pairs = [("R3", p) for p in sorted(A.participant_id)[:3]]
    pairs += [("R10", s.replace("sub-", "")) for s in rep["fcp_R10"]["first_subjects"][:2]]
    eq = []
    for rel, pid in pairs:
        a = sha256_url(f"{FCP}/{FCP_PREFIX}/cmi_bids_{rel}/{set_path(pid)}")
        b = sha256_url(f"{ON}/{DS[rel]}/{set_path(pid)}")
        eq.append(dict(release=rel, participant_id=pid, fcp=a, openneuro=b,
                       identical=(a.get("sha256") is not None and a.get("sha256") == b.get("sha256"))))
    rep["equivalence"] = eq
    os.makedirs(OUT, exist_ok=False)
    with open(f"{OUT}/diagnostic.json", "w") as fh:
        json.dump(rep, fh, indent=1)
    # printed summary, read back from the file
    with open(f"{OUT}/diagnostic.json") as fh:
        R = json.load(fh)
    print(f"output: {OUT}/diagnostic.json\n")
    for rel in ("R10", "R11"):
        f = R[f"fcp_{rel}"]
        print(f"FCP-INDI cmi_bids_{rel}: listing HTTP {f['listing_status']}, {f['n_subject_folders']} subject folders "
              f"(first: {f['first_subjects']}); root files: {f['root_files']}; other folders: {f['other_folders']}")
        o = R[f"openneuro_{rel}_participants"]
        print(f"OpenNeuro {DS[rel]} participants.tsv: {o}")
    print("\nEQUIVALENCE (FCP-INDI vs OpenNeuro, resting-state .set):")
    for e in R["equivalence"]:
        print(f"  {e['release']} {e['participant_id']}: identical={e['identical']} | "
              f"FCP {e['fcp'].get('status')} {e['fcp'].get('bytes')} B {str(e['fcp'].get('sha256'))[:12]} | "
              f"OpenNeuro {e['openneuro'].get('status')} {e['openneuro'].get('bytes')} B {str(e['openneuro'].get('sha256'))[:12]}")


if __name__ == "__main__":
    main()
