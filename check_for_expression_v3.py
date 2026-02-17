#!/usr/bin/env python3

import pandas as pd
import glob
import sys

# ----------------------------
# Parse command line arguments
# ----------------------------
if len(sys.argv) != 4:
    sys.exit(
        "Usage:\n"
        "  python script.py <pangenome_csv> <tsv_glob_pattern> <output_csv>\n\n"
        "Example:\n"
        "  python script.py gene_presence_absence.csv '/path/to/*301.tsv' out.csv"
    )

pangenome_csv = sys.argv[1]
tsv_pattern = sys.argv[2]
output_csv = sys.argv[3]


# ----------------------------
# Load pangenome table
# ----------------------------
df = pd.read_csv(pangenome_csv, dtype=str, low_memory=False).fillna("")


# ----------------------------
# Collect expressed labels
# ----------------------------
expressed_labels = set()

tsv_files = glob.glob(tsv_pattern)
if not tsv_files:
    sys.exit(f"ERROR: No TSV files matched pattern: {tsv_pattern}")

for tsv_file in tsv_files:
    with open(tsv_file, "r", encoding="utf-8") as f:
        first = True
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if not cols:
                continue

            label = cols[0].strip()
            if not label:
                continue

            # skip header row
            if first and label.lower() in {
                "locus_tag", "locus", "label", "id", "protein", "protein_id"
            }:
                first = False
                continue

            first = False
            expressed_labels.add(label)


# ----------------------------
# Expression filter function
# ----------------------------
def check_expression_cell(cell):
    cell = (cell or "").strip()
    if cell == "":
        return ""

    parts = [lbl.strip() for lbl in cell.split(";") if lbl.strip()]
    kept = [
        lbl for lbl in parts
        if "refound" not in lbl.lower() and lbl in expressed_labels
    ]

    return ";".join(kept) if kept else "NA"


# ----------------------------
# Apply to isolate columns
# ----------------------------
meta_cols = ["Gene", "Non-unique Gene name", "Annotation"]
isolate_cols = [c for c in df.columns if c not in meta_cols]

df_checked = df.copy()

for c in isolate_cols:
    df_checked[c] = df_checked[c].apply(check_expression_cell)


# ----------------------------
# Write output
# ----------------------------
df_checked.to_csv(output_csv, index=False)

print(f"Wrote: {output_csv}")
print(f"Loaded TSV files: {len(tsv_files)}")
print(f"Unique expressed labels: {len(expressed_labels)}")
