#!/usr/bin/env python3
"""
Build a cCRE × fHSS × HAR × HAQER table 
for fixed HAR and HAQER regions.
"""

import subprocess
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent

CCRE_TSV  = BASE / "ccre.9num1kb.lss.header.tsv"
FHSS_BED  = BASE / "9way.human.hprcv2rename.allchr.bed"
HAR_BED   = BASE / "Fixed_HAR_merged.bed"
HAQER_BED = BASE / "Fixed_HAQER2025_merged.bed"

OUTDIR = BASE / "ccre_har_haqer"
OUTCSV = OUTDIR / "ccre_har_haqer_fix.csv"


FINAL_COLS = [
    "chr", "start", "end", "rDHS", "cCRE_ID", "cCRE_class", "Gene_IDs",
    "Human", "Hominini", "Homininae", "Hominidae", "Hominoidea", "Conserved",
    "ccre_length", "fhss_count", "pct_fhss", "fhss_category",
    "har_count", "pct_har", "haqer_count", "pct_haqer",
]


def run_bedtools_count(ccre_bed_path, feature_bed_path):
    """bedtools intersect -a ccre.bed -b feature.bed -c → dict[idx] = count."""
    cmd = ["bedtools", "intersect",
           "-a", str(ccre_bed_path),
           "-b", str(feature_bed_path),
           "-c"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    counts = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        idx = int(fields[3])
        counts[idx] = int(fields[-1])
    return counts


def run_bedtools_overlap(ccre_bed_path, feature_bed_path):
    """bedtools intersect -wao → dict[idx] = total overlap bp."""
    cmd = ["bedtools", "intersect",
           "-a", str(ccre_bed_path),
           "-b", str(feature_bed_path),
           "-wao"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    overlaps = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        idx = int(fields[3])
        bp = int(fields[-1])
        overlaps[idx] = overlaps.get(idx, 0) + bp
    return overlaps


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("Reading cCRE table...")
    df = pd.read_csv(CCRE_TSV, sep="\t", encoding="utf-8-sig")
    n = len(df)
    print(f"  {n:,} cCREs loaded")

    df["ccre_length"] = df["end"] - df["start"]

    # Temp BED with row index in column 4 for stable lookup
    tmp_dir = BASE / "tmp_ccre_har_haqer"
    tmp_dir.mkdir(exist_ok=True)
    tmp_bed = tmp_dir / "ccre_temp.bed"
    with open(tmp_bed, "w") as fh:
        for i, row in enumerate(df[["chr", "start", "end"]].itertuples(index=False)):
            fh.write(f"{row.chr}\t{row.start}\t{row.end}\t{i}\n")

    # ── fHSS counts ────────────────────────────────────────────────────────
    print("Counting fHSS variants per cCRE (bedtools intersect -c)...")
    fhss_counts = run_bedtools_count(tmp_bed, FHSS_BED)
    df["fhss_count"] = [fhss_counts.get(i, 0) for i in range(n)]
    print(f"  {(df['fhss_count'] > 0).sum():,} cCREs have ≥1 fHSS")

    df["pct_fhss"] = 0.0
    mask = df["ccre_length"] > 0
    df.loc[mask, "pct_fhss"] = df.loc[mask, "fhss_count"] / df.loc[mask, "ccre_length"] * 100

    df["fhss_category"] = pd.cut(
        df["pct_fhss"],
        bins=[-1, 0, 10, 25, 50, 75, 101],
        labels=["0%", "1-10%", "11-25%", "26-50%", "51-75%", ">75%"],
    )

    # ── HAR overlaps ──────────────────────────────────────────────────────
    print("Counting HAR overlap bp per cCRE (bedtools intersect -wao)...")
    har_overlaps = run_bedtools_overlap(tmp_bed, HAR_BED)
    df["har_count"] = [har_overlaps.get(i, 0) for i in range(n)]
    df["pct_har"] = (df["har_count"] / df["ccre_length"]) * 100
    print(f"  {(df['har_count'] > 0).sum():,} cCREs overlap HARs")

    # ── HAQER overlaps ────────────────────────────────────────────────────
    print("Counting HAQER overlap bp per cCRE (bedtools intersect -wao)...")
    haqer_overlaps = run_bedtools_overlap(tmp_bed, HAQER_BED)
    df["haqer_count"] = [haqer_overlaps.get(i, 0) for i in range(n)]
    df["pct_haqer"] = (df["haqer_count"] / df["ccre_length"]) * 100
    print(f"  {(df['haqer_count'] > 0).sum():,} cCREs overlap HAQERs")

    tmp_bed.unlink(missing_ok=True)
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    out = df[FINAL_COLS]
    print(f"Writing {OUTCSV} ...")
    out.to_csv(OUTCSV, index=False)
    print(f"  {len(out):,} rows × {len(out.columns)} cols written")


if __name__ == "__main__":
    main()
