#!/usr/bin/env python3
"""
prepare_btcr_links.py
Create standardized FASTQ symlinks for bulk BCR/TCR samples.
Reads samples.csv and creates bcl2fastq-style naming in output directory.
"""
import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare BTCRFlow FASTQ links")
    parser.add_argument("--sample-csv", required=True, help="Path to samples.csv")
    parser.add_argument("--output-dir", required=True, help="Output link directory")
    parser.add_argument("--manifest", required=True, help="Output manifest JSON path")
    return parser.parse_args()


def find_fastq(raw_path, index, read_num):
    """Find FASTQ file for given read number."""
    patterns = [
        f"*{index}*R{read_num}*.fastq.gz",
        f"*{index}*R{read_num}*.fq.gz",
        f"*{index}*_{read_num}.fastq.gz",
        f"*{index}*_{read_num}.fq.gz",
    ]
    for pat in patterns:
        matches = sorted(glob.glob(os.path.join(raw_path, pat)))
        if matches:
            return matches[0]
    return None


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    manifest = {"samples": {}}

    with open(args.sample_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["sample"]
            index = row.get("index", sample)
            raw_path = row["path"]
            platform = row.get("platform", "tcr").lower()

            sample_dir = os.path.join(args.output_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)

            r1_src = find_fastq(raw_path, index, 1)
            r2_src = find_fastq(raw_path, index, 2)

            if not r1_src or not r2_src:
                print(f"WARNING: Missing FASTQ for {sample} in {raw_path}", file=sys.stderr)
                continue

            r1_link = os.path.join(sample_dir, f"{sample}_S1_L001_R1_001.fastq.gz")
            r2_link = os.path.join(sample_dir, f"{sample}_S1_L001_R2_001.fastq.gz")

            if not os.path.exists(r1_link):
                os.symlink(os.path.abspath(r1_src), r1_link)
            if not os.path.exists(r2_link):
                os.symlink(os.path.abspath(r2_src), r2_link)

            manifest["samples"][sample] = {
                "platform": platform,
                "r1": os.path.abspath(r1_link),
                "r2": os.path.abspath(r2_link),
                "raw_path": raw_path,
                "index": index,
            }

    with open(args.manifest, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest written: {args.manifest}")
    print(f"Samples linked: {len(manifest['samples'])}")


if __name__ == "__main__":
    main()
