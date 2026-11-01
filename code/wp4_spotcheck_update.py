# WP4 step 3 of 3: apply the author's spot-check decisions (26 Oct 2026) to audit_table.csv.
# Run order: wp4_make_audit_tables.py -> wp4_a3_update.py -> wp4_spotcheck_update.py (same folder).
import pandas as pd
t = pd.read_csv("audit_table.csv")
t.loc[t.audit_id == "A05", "order_category"] = "Not determinable"
t.loc[t.audit_id == "A05", "notes"] = ("author spot-check 26 Oct 2026: sentence order unclear -> Not determinable "
                                       "(was References before detecting)")
t.loc[t.audit_id == "A04", "notes"] = t.loc[t.audit_id == "A04", "notes"].astype(str) + "; author spot-check: classification confirmed"
for i in ("A10", "A12"):
    t.loc[t.audit_id == i, "order_category"] = "Not accessible"
    t.loc[t.audit_id == i, "notes"] = "full text not accessible to Claude or author"
t.to_csv("audit_table.csv", index=False)
print(t.order_category.value_counts())
