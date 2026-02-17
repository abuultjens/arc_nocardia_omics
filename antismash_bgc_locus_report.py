#!/usr/bin/env python3
"""
antiSMASH BGC locus-tag prevalence + BiG-SCAPE GCF mapping (all classes) + Panaroo orthogroup linking + matrices,
restricted to a user-specified isolate subset.

NEW:
- Two positional args:
  1) gene_presence_absence.csv (Panaroo output)
  2) isolate_file_list.txt      (one filename per line, e.g. SP0015; no header)

The script will only process isolates listed in isolate_file_list.txt.
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

try:
    import pandas as pd
except ImportError:
    raise SystemExit("ERROR: This script requires pandas. Install with: pip install pandas")


# ----------------------------
# Helpers: isolate list parsing
# ----------------------------

def load_isolate_ids_from_file(list_path: Path) -> List[str]:
    """
    Read a 1-column text file containing filenames (e.g. SP0015.gff) or isolate IDs (e.g. SP0015),
    one per line, no header. Returns a de-duplicated list of isolate IDs preserving order.
    """
    if not list_path.exists():
        raise SystemExit(f"ERROR: isolate file list not found: {list_path}")

    ids: List[str] = []
    seen = set()

    with list_path.open("r", errors="replace") as fh:
        for raw in fh:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                continue
            # Take basename only in case paths are present
            s = s.replace("\\", "/").split("/")[-1]
            # If it's a filename, strip extension(s)
            s = re.sub(r"\.(gff|gff3|gbk|gb|fna|faa|fa|fasta|tsv|csv|txt)$", "", s, flags=re.IGNORECASE)
            # Now s should be isolate ID like SP0015
            if s not in seen:
                ids.append(s)
                seen.add(s)

    if not ids:
        raise SystemExit(f"ERROR: No isolate IDs found in {list_path}")

    return ids


# ----------------------------
# antiSMASH GBK parsing helpers
# ----------------------------

_LOCUS_TAG_RE = re.compile(r'^\s*/locus_tag="([^"]+)"\s*$')


def extract_locus_tags_from_gbk(gbk_path: Path) -> Set[str]:
    """Lightweight extraction of locus_tag qualifiers from a GenBank file."""
    loci = set()
    with gbk_path.open("r", errors="replace") as fh:
        for line in fh:
            m = _LOCUS_TAG_RE.match(line.rstrip("\n"))
            if m:
                loci.add(m.group(1))
    return loci


def parse_region_location_and_products(gbk_path: Path) -> Tuple[Optional[Tuple[int, int]], List[str]]:
    """
    Parse region coordinates and product types from an antiSMASH region GenBank file.

    Returns:
      (start, end) tuple (1-based inclusive) if detected else None,
      list of product types (unique, order preserved) if detected else [].
    """
    products: List[str] = []
    seen_prod = set()

    region_loc: Optional[Tuple[int, int]] = None

    in_features = False
    current_feature = None
    current_loc = None

    feat_line_re = re.compile(r"^\s{5}(\S+)\s+(.+)\s*$")
    prod_re = re.compile(r'^\s*/product="([^"]+)"\s*$')

    def loc_to_span(loc_str: str) -> Optional[Tuple[int, int]]:
        nums = [int(x) for x in re.findall(r"(\d+)", loc_str)]
        if len(nums) >= 2:
            return (min(nums), max(nums))
        return None

    with gbk_path.open("r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            if line.startswith("ORIGIN"):
                break

            m = feat_line_re.match(line)
            if m:
                current_feature = m.group(1)
                current_loc = m.group(2).strip()

                if current_feature.lower() == "region":
                    span = loc_to_span(current_loc)
                    if span:
                        region_loc = span
                continue

            pm = prod_re.match(line.strip())
            if pm and current_feature is not None:
                if current_feature.lower() in ("region", "protocluster"):
                    prod = pm.group(1).strip()
                    if prod and prod not in seen_prod:
                        products.append(prod)
                        seen_prod.add(prod)

    if region_loc is None:
        # fallback: protocluster span
        with gbk_path.open("r", errors="replace") as fh:
            in_features = False
            current_feature = None
            current_loc = None
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("FEATURES"):
                    in_features = True
                    continue
                if not in_features:
                    continue
                if line.startswith("ORIGIN"):
                    break

                m = feat_line_re.match(line)
                if m:
                    current_feature = m.group(1)
                    current_loc = m.group(2).strip()
                    if current_feature.lower() == "protocluster":
                        span = loc_to_span(current_loc)
                        if span:
                            region_loc = span
                            break

    return region_loc, products


def region_size_kb_from_span(span: Optional[Tuple[int, int]]) -> Optional[float]:
    if span is None:
        return None
    start, end = span
    if start <= 0 or end <= 0 or end < start:
        return None
    bp = end - start + 1
    return round(bp / 1000.0, 3)


# ----------------------------
# BiG-SCAPE clustering parsing (all classes)
# ----------------------------

def norm_bgc_id(x: str) -> str:
    """
    Normalise BGC IDs to match antiSMASH region basenames:
      SP0040_contig_1.region007.gbk -> SP0040_contig_1.region007
      /path/to/SP0040_contig_1.region007.gbk -> SP0040_contig_1.region007
    """
    x = str(x).strip()
    x = x.replace("\\", "/")
    x = x.split("/")[-1]
    for ext in [".gbk", ".gbff", ".genbank", ".gb", ".json"]:
        if x.endswith(ext):
            x = x[: -len(ext)]
    return x


def parse_bigscape_clustering_file(path: Path) -> List[Tuple[str, str]]:
    """
    Parse BiG-SCAPE *_clustering_c*.tsv.

    Returns list of (bgc_name, family_number) strings.
    Handles whitespace-delimited content robustly.
    """
    pairs: List[Tuple[str, str]] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            bgc = parts[0]
            fam = parts[1]
            pairs.append((bgc, fam))
    return pairs


def load_bigscape_all_classes(run_dir: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """
    Load mapping from ALL class clustering tables under a BiG-SCAPE run directory.

    Returns mapping:
      region_id -> {
        "bigscape_class": <class>,
        "family_number": <family>,
        "gcf_id": <class>_GCF<family>
      }
    """
    mapping: Dict[str, Dict[str, str]] = {}

    if run_dir is None:
        return mapping
    if not run_dir.exists():
        print(f"WARNING: BiG-SCAPE run directory not found: {run_dir} (GCF IDs will be blank)")
        return mapping

    clustering_files = sorted(run_dir.glob("*/*_clustering_c*.tsv"))
    if not clustering_files:
        clustering_files = sorted(run_dir.glob("**/*_clustering_c*.tsv"))

    if not clustering_files:
        print(f"WARNING: No *_clustering_c*.tsv files found under {run_dir} (GCF IDs will be blank)")
        return mapping

    n_pairs = 0
    for f in clustering_files:
        bigscape_class = f.parent.name
        pairs = parse_bigscape_clustering_file(f)
        n_pairs += len(pairs)

        for bgc_raw, fam_raw in pairs:
            bgc = norm_bgc_id(bgc_raw)
            fam = str(fam_raw).strip()
            if not bgc or not fam:
                continue
            gcf_id = f"{bigscape_class}_GCF{fam}"
            mapping[bgc] = {
                "bigscape_class": bigscape_class,
                "family_number": fam,
                "gcf_id": gcf_id
            }

    print(f"Loaded BiG-SCAPE clustering from: {run_dir}")
    print(f"Found {len(clustering_files)} clustering tables, parsed {n_pairs} rows, mapped {len(mapping)} BGC IDs.")
    return mapping


# ----------------------------
# Panaroo parsing helpers
# ----------------------------

def load_panaroo_locus_to_group(gpa_csv: Path, allowed_isolates: Set[str]) -> Dict[str, str]:
    """
    Load Panaroo gene_presence_absence.csv and build mapping:
      locus_tag -> panaroo_gene_cluster_id  (value from the 'Gene' column)

    Only uses isolate columns that match allowed_isolates (intersection).
    Splits multi-locus cells on ';' or ',' only.
    """
    if not gpa_csv.exists():
        raise SystemExit(f"ERROR: gene_presence_absence.csv not found: {gpa_csv}")

    df = pd.read_csv(gpa_csv, dtype=str).fillna("")

    if "Gene" not in df.columns:
        raise SystemExit("ERROR: gene_presence_absence.csv does not contain a 'Gene' column.")

    fixed_cols = {"Gene", "Non-unique Gene name", "Annotation"}
    all_isolate_cols = [c for c in df.columns if c not in fixed_cols]

    # Restrict to requested isolate subset
    isolate_cols = [c for c in all_isolate_cols if c in allowed_isolates]

    missing = sorted(list(allowed_isolates - set(all_isolate_cols)))
    if missing:
        print(f"WARNING: {len(missing)} requested isolates were not found as columns in gene_presence_absence.csv.")
        print("         Example missing:", ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))

    if not isolate_cols:
        raise SystemExit("ERROR: None of the requested isolates were found as columns in gene_presence_absence.csv.")

    locus_to_group: Dict[str, str] = {}
    dup_loci = 0

    for _, row in df.iterrows():
        group_id = str(row.get("Gene", "")).strip()
        if not group_id:
            continue

        for col in isolate_cols:
            cell = str(row.get(col, "")).strip()
            if not cell:
                continue

            parts = re.split(r"[;,]\s*", cell)
            parts = [p.strip() for p in parts if p.strip()]

            for locus in parts:
                if locus in locus_to_group and locus_to_group[locus] != group_id:
                    dup_loci += 1
                    continue
                locus_to_group.setdefault(locus, group_id)

    if dup_loci > 0:
        print(f"WARNING: {dup_loci} locus tags were found assigned to multiple Panaroo groups. Kept first assignment.")

    print(f"Loaded Panaroo groups from: {gpa_csv}")
    print(f"Using {len(isolate_cols)} isolate columns from Panaroo table (restricted subset).")
    print(f"Mapped {len(locus_to_group)} locus tags to Panaroo gene clusters.")
    return locus_to_group


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene_presence_absence_csv",
                    help="Panaroo gene_presence_absence.csv (maps locus tags to orthogroup IDs)")
    ap.add_argument("isolate_file_list",
                    help="1-col text file (no header) listing isolate filenames (e.g. SP0015.gff) or isolate IDs (SP0015)")
    ap.add_argument("--antismash_dir", required=True,
                    help="Parent directory containing SP**** antiSMASH folders")
    ap.add_argument("--bigscape_run_dir", default=None,
                    help="BiG-SCAPE run directory containing class folders (e.g. .../<timestamp>_hybrids_glocal/)")
    ap.add_argument("--out_prefix", default="antismash_bgc",
                    help="Output prefix for CSV files")
    args = ap.parse_args()

    gpa_csv = Path(args.gene_presence_absence_csv)
    isolate_list_path = Path(args.isolate_file_list)

    antismash_base = Path(args.antismash_dir)
    bigscape_run_dir = Path(args.bigscape_run_dir) if args.bigscape_run_dir else None

    if not antismash_base.exists():
        raise SystemExit(f"ERROR: antismash_dir not found: {antismash_base}")

    isolate_ids = load_isolate_ids_from_file(isolate_list_path)
    isolate_set = set(isolate_ids)

    # Restrict isolate dirs to requested subset
    isolate_dirs = []
    for iso in isolate_ids:
        d = antismash_base / iso
        if d.is_dir():
            isolate_dirs.append(d)
        else:
            print(f"WARNING: isolate directory not found under --antismash_dir: {d}")

    if not isolate_dirs:
        raise SystemExit("ERROR: None of the requested isolate directories were found under --antismash_dir.")

    print(f"Requested isolates: {len(isolate_ids)}")
    print(f"Found isolate folders: {len(isolate_dirs)} under {antismash_base}")

    # Load mappings
    locus_to_panaroo_group = load_panaroo_locus_to_group(gpa_csv, allowed_isolates=isolate_set)
    bgc_to_gcf_info = load_bigscape_all_classes(bigscape_run_dir)

    detail_rows: List[Dict[str, object]] = []

    locus_to_isolates: Dict[str, Set[str]] = defaultdict(set)
    locus_to_regions: Dict[str, Set[str]] = defaultdict(set)
    locus_to_products: Dict[str, Set[str]] = defaultdict(set)
    locus_to_gcf: Dict[str, Set[str]] = defaultdict(set)
    locus_to_sizes_kb: Dict[str, List[float]] = defaultdict(list)
    locus_to_bigscape_classes: Dict[str, Set[str]] = defaultdict(set)
    locus_to_family_numbers: Dict[str, Set[str]] = defaultdict(set)
    locus_to_panaroo_groups: Dict[str, Set[str]] = defaultdict(set)

    locus_seen_anywhere: Dict[str, Set[str]] = defaultdict(set)
    locus_seen_in_bgc: Dict[str, Set[str]] = defaultdict(set)

    isolate_to_bgc_loci: Dict[str, Set[str]] = defaultdict(set)

    for iso_dir in isolate_dirs:
        isolate = iso_dir.name

        region_files = sorted(iso_dir.glob("*region*.gbk"))
        if not region_files:
            print(f"WARNING: No *region*.gbk files found in {iso_dir}")
            continue

        for region_gbk in region_files:
            region_id = region_gbk.name.replace(".gbk", "")
            span, products = parse_region_location_and_products(region_gbk)
            size_kb = region_size_kb_from_span(span)

            gcf_info = bgc_to_gcf_info.get(region_id, {})
            bigscape_class = gcf_info.get("bigscape_class", "")
            family_number = gcf_info.get("family_number", "")
            gcf_id = gcf_info.get("gcf_id", "")

            loci = extract_locus_tags_from_gbk(region_gbk)
            for locus in loci:
                isolate_to_bgc_loci[isolate].add(locus)

                locus_seen_in_bgc[locus].add(isolate)
                locus_seen_anywhere[locus].add(isolate)

                locus_to_isolates[locus].add(isolate)
                locus_to_regions[locus].add(f"{isolate}:{region_id}")

                for p in products:
                    locus_to_products[locus].add(p)

                if gcf_id:
                    locus_to_gcf[locus].add(gcf_id)
                if bigscape_class:
                    locus_to_bigscape_classes[locus].add(bigscape_class)
                if family_number:
                    locus_to_family_numbers[locus].add(family_number)

                if size_kb is not None:
                    locus_to_sizes_kb[locus].append(size_kb)

                panaroo_group = locus_to_panaroo_group.get(locus, "")
                if panaroo_group:
                    locus_to_panaroo_groups[locus].add(panaroo_group)

                detail_rows.append({
                    "isolate": isolate,
                    "region_id": region_id,
                    "locus_tag": locus,
                    "panaroo_gene_cluster": panaroo_group,
                    "antismash_product_types": ";".join(products) if products else "",
                    "region_size_kb": size_kb if size_kb is not None else "",
                    "bigscape_class": bigscape_class,
                    "bigscape_family_number": family_number,
                    "bigscape_gcf_id": gcf_id
                })

        # Optional: scan main genome GBK to detect loci outside BGCs
        main_gbk = iso_dir / f"{isolate}.gbk"
        if main_gbk.exists():
            all_loci = extract_locus_tags_from_gbk(main_gbk)
            for locus in all_loci:
                locus_seen_anywhere[locus].add(isolate)

    # ---------------------------------
    # Write isolate-region-locus mapping
    # ---------------------------------
    detail_df = pd.DataFrame(detail_rows)
    detail_out = f"{args.out_prefix}_isolate_region_locus.csv"
    detail_df.to_csv(detail_out, index=False)

    # ----------------------------
    # Build locus summary dataframe
    # ----------------------------
    summary_rows = []
    for locus in sorted(locus_to_isolates.keys()):
        isolates = sorted(locus_to_isolates[locus])
        regions = sorted(locus_to_regions[locus])
        products = sorted(locus_to_products.get(locus, set()))
        gcfs = sorted(locus_to_gcf.get(locus, set()))
        bs_classes = sorted(locus_to_bigscape_classes.get(locus, set()))
        fams = sorted(locus_to_family_numbers.get(locus, set()))
        pan_groups = sorted(locus_to_panaroo_groups.get(locus, set()))
        sizes = locus_to_sizes_kb.get(locus, [])

        size_min = min(sizes) if sizes else ""
        size_median = (pd.Series(sizes).median() if sizes else "")
        size_max = max(sizes) if sizes else ""

        ever_outside = (locus_seen_anywhere[locus] != locus_seen_in_bgc[locus])

        summary_rows.append({
            "locus_tag": locus,
            "panaroo_gene_clusters": ";".join(pan_groups),
            "num_isolates_in_bgc": len(isolates),
            "isolates_in_bgc": ";".join(isolates),
            "regions_seen_in": ";".join(regions),
            "antismash_product_types": ";".join(products),
            "region_size_kb_min": size_min,
            "region_size_kb_median": size_median,
            "region_size_kb_max": size_max,
            "bigscape_classes": ";".join(bs_classes),
            "bigscape_family_numbers": ";".join(fams),
            "bigscape_gcf_ids": ";".join(gcfs),
            "ever_seen_outside_bgc": ever_outside
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_out = f"{args.out_prefix}_locus_summary.csv"
    summary_df.to_csv(summary_out, index=False)

    # -----------------------------------------
    # Presence/absence matrix (isolates ?? loci)
    # -----------------------------------------
    isolates_all = [d.name for d in isolate_dirs]  # preserve requested subset order (minus missing)
    loci_all = sorted(locus_to_isolates.keys())

    matrix = []
    for iso in isolates_all:
        present = isolate_to_bgc_loci.get(iso, set())
        row = {"isolate": iso}
        for locus in loci_all:
            row[locus] = 1 if locus in present else 0
        matrix.append(row)

    matrix_df = pd.DataFrame(matrix)
    matrix_out = f"{args.out_prefix}_presence_absence_matrix.csv"
    matrix_df.to_csv(matrix_out, index=False)

    # -----------------------
    # Heatmap-ready long table
    # -----------------------
    long_df = matrix_df.melt(id_vars=["isolate"], var_name="locus_tag", value_name="present")

    meta_cols = [
        "locus_tag",
        "panaroo_gene_clusters",
        "antismash_product_types",
        "bigscape_gcf_ids",
        "bigscape_classes",
        "num_isolates_in_bgc",
        "region_size_kb_median",
    ]
    meta_df = summary_df[meta_cols].copy()
    heatmap_df = long_df.merge(meta_df, on="locus_tag", how="left")

    heatmap_out = f"{args.out_prefix}_heatmap_long.csv"
    heatmap_df.to_csv(heatmap_out, index=False)

    print("\nDone!")
    print(f"Wrote: {detail_out}")
    print(f"Wrote: {summary_out}")
    print(f"Wrote: {matrix_out}")
    print(f"Wrote: {heatmap_out}")


if __name__ == "__main__":
    main()
