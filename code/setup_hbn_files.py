# =============================================================================
# SETUP: make sure the three HBN R3 RestingState recordings are present + verified
# Run this cell FIRST in every new Colab session (Colab wipes /content on restart).
#
# 1. Searches for existing copies anywhere under /content (incl. a mounted Drive).
# 2. Otherwise downloads from the public FCP-INDI S3 bucket over plain HTTPS
#    (uses `requests`, preinstalled in Colab -- no s3fs / awscli install needed).
# 3. Verifies each .set against the SHA-256 recorded in the handover.
# Never overwrites a file that already passes verification.
# =============================================================================
import os, glob, shutil, hashlib, requests

S3_HTTPS = "https://fcp-indi.s3.amazonaws.com/data/Projects/HBN/BIDS_EEG/cmi_bids_R3"
PARTICIPANTS = {   # id : (target folder, expected SHA-256 of the .set, from the handover)
    "NDARAA948VFH": ("/content/hbn_one",
                     "e7f078a9d5406589d904d29a46c784aab5d8ed1fb09c9a9ad7a71dbebd02f121"),
    "NDARAC853DTE": ("/content/hbn_two",
                     "a43f1c48ec64c0799d8cf6853ff5e4c0970aa7d9049b8ca337430d4b05733cad"),
    "NDARAD774HAZ": ("/content/hbn_three",
                     "ccae811acbc740709d348b6ffb3c0876c8f8872a4de77c4ce8f1d09520a542cd"),
}
SUFFIXES = ["_eeg.set", "_events.tsv", "_channels.tsv", "_eeg.json"]


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def find_existing(fname):
    """Search /content recursively (covers a mounted Google Drive at /content/drive)."""
    hits = glob.glob(f"/content/**/{fname}", recursive=True)
    return hits[0] if hits else None


def download(url, dest):
    tmp = dest + ".part"
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        with open(tmp, "wb") as f:
            for b in r.iter_content(chunk_size=1 << 20):
                f.write(b)
                done += len(b)
        # content-length is the size ON THE WIRE; if the server compressed the response,
        # requests decompresses it, so sizes legitimately differ. Compare only when the
        # response is not content-encoded. The SHA-256 check is the authoritative test.
        encoded = r.headers.get("content-encoding", "").lower() not in ("", "identity")
        if total and not encoded and done != total:
            raise IOError(f"incomplete download: {done} of {total} bytes")
    os.replace(tmp, dest)                      # only replace once fully written


report = []
for pid, (folder, expected) in PARTICIPANTS.items():
    os.makedirs(folder, exist_ok=True)
    for sfx in SUFFIXES:
        fname = f"sub-{pid}_task-RestingState{sfx}"
        dest = os.path.join(folder, fname)
        if os.path.exists(dest):
            source = "already present"
        else:
            found = find_existing(fname)
            if found and os.path.abspath(found) != os.path.abspath(dest):
                shutil.copy2(found, dest)
                source = f"copied from {found}"
            else:
                download(f"{S3_HTTPS}/sub-{pid}/eeg/{fname}", dest)
                source = "downloaded from S3"
        size = os.path.getsize(dest)
        status = ""
        if sfx == "_eeg.set":
            got = sha256(dest)
            status = "SHA-256 OK" if got == expected else f"SHA-256 MISMATCH (got {got[:16]}...)"
        report.append((pid, fname, size, source, status))

print(f"{'participant':<14}{'file':<52}{'bytes':>12}  source / check")
for pid, fname, size, source, status in report:
    print(f"{pid:<14}{fname:<52}{size:>12,}  {source}  {status}")

bad = [r for r in report if "MISMATCH" in r[4]]
if bad:
    raise RuntimeError("Checksum mismatch -- do NOT run the diagnostic until resolved.")
print("\nAll three recordings present and verified. Now run the QC diagnostic cell.")
