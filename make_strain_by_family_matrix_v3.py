#!/usr/bin/env python3

import sys
import re
import pandas as pd
from pathlib import Path


def sanitise_type_label(s):
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s)
    s = s.strip("_")
    return s if s else "UNKNOWN"


def normalise_family_number(x):
    if pd.isna(x):
        return None
    try:
        fx = float(x)
        if fx.is_integer():
            return str(int(fx))
        return str(x).strip()
    except Exception:
        return str(x).strip()


def parse_clustering_file(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    if "#BGC Name" not in df.columns or "Family Number" not in df.columns:
        raise ValueError(
            "Expected '#BGC Name' and 'Family Number' in {}. Got {}".format(
                path, df.columns.tolist()
            )
        )

    type_label = sanitise_type_label(path.parent.name)

    df = df[df["#BGC Name"].astype(str).str.startswith("SP")].copy()

    df["Strain"] = df["#BGC Name"].astype(str).str.split("_").str[0]

    df["FamilyNumberNorm"] = df["Family Number"].apply(normalise_family_number)

    df["Feature"] = type_label + "_" + df["FamilyNumberNorm"].astype(str)

    df = df[["Strain", "Feature"]].dropna().drop_duplicates()

    return df


def main(out_path, files):
    all_mappings = []

    for f in files:
        try:
            p = Path(f)
            parsed = parse_clustering_file(p)
            all_mappings.append(parsed)
            print("[OK] Parsed {} (type={})".format(
                f, sanitise_type_label(p.parent.name)))
        except Exception as e:
            print("[WARN] Skipping {}: {}".format(f, e))

    if not all_mappings:
        sys.exit("[ERROR] No valid clustering files parsed!")

    merged = pd.concat(all_mappings, ignore_index=True)

    print("[INFO] Example features:",
          merged["Feature"].dropna().unique()[:10])

    merged["value"] = 1

    matrix = merged.pivot_table(
        index="Strain",
        columns="Feature",
        values="value",
        aggfunc="max",
        fill_value=0
    )

    matrix.columns = matrix.columns.astype(str)

    matrix.sort_index(axis=0, inplace=True)
    matrix.sort_index(axis=1, inplace=True)

    matrix.to_csv(out_path)

    print("[OK] Wrote matrix:", out_path)
    print("[INFO] Shape: strains={}, features={}".format(
        matrix.shape[0], matrix.shape[1]))


if __name__ == "__main__":

    if len(sys.argv) < 3:
        sys.exit(
            "Usage:\n"
            "  python make_strain_by_family_matrix.py OUTPUT.csv */*clustering_c0.30.tsv"
        )

    out_file = sys.argv[1]
    clustering_files = sys.argv[2:]

    main(out_file, clustering_files)
