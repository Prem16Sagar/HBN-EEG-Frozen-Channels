# Builds the Release 10 census tables (recordings, channels, episodes) from the per-child
# JSON files written by a2_r10r11_census.py, using the table-building code of wp1_census.py
# itself (copied at run time, not rewritten), so the format equals the Release 1-9 tables.
# Usage: put wp1_census.py beside this script and the per_child folder
# (revision_B/census/census_v1_r10r11/per_child) in the working directory, then run it.
import glob, json, os, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "wp1_census.py")).read()
a = src.index("    ch, ep, rec = [], [], []"); b = src.index("    CH, EP, REC = pd.DataFrame(ch)")
J = [json.load(open(f)) for f in sorted(glob.glob("per_child/R10_*.json")) if not f.endswith((".nofile.json", ".failed.json"))]
ns = {"J": J, "KNOWN": 193893.78125}
exec("\n".join(line[4:] for line in src[a:b].splitlines()), ns)
for name, rows in (("recordings", ns["rec"]), ("channels", ns["ch"]), ("episodes", ns["ep"])):
    pd.DataFrame(rows).to_csv(f"census_{name}_R10.csv", index=False)
print(f"{len(J)} Release 10 recordings written")
