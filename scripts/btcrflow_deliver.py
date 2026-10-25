#!/usr/bin/env python3
"""
btcrflow_deliver.py
Generate delivery manifest and MD5 checksums for BTCRFlow outputs.
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="BTCRFlow delivery manifest generator")
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--output-dir", required=True, help="data_deliver directory")
    parser.add_argument("--manifest", required=True, help="Output manifest JSON path")
    parser.add_argument("--md5", required=True, help="Output MD5 file path")
    return parser.parse_args()


def md5_file(filepath, chunk_size=8192):
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def collect_delivery_files(base_dir):
    """Collect key output files for delivery."""
    files = []
    patterns = [
        ("01.qc/multiqc/multiqc.html", "QC Report"),
        ("03.mixcr/postanalysis.json.gz", "MiXCR Postanalysis"),
    ]
    for pattern, desc in patterns:
        full = os.path.join(base_dir, pattern)
        if os.path.exists(full):
            files.append({"path": pattern, "description": desc})

    mixcr_dir = os.path.join(base_dir, "03.mixcr")
    if os.path.isdir(mixcr_dir):
        for f in sorted(os.listdir(mixcr_dir)):
            if f.endswith(".clns"):
                files.append({
                    "path": f"03.mixcr/{f}",
                    "description": f"MiXCR clones: {f}",
                })

    return files


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.manifest), exist_ok=True)

    workdir = os.getcwd()
    files = collect_delivery_files(workdir)

    manifest = {
        "project_name": args.project_name,
        "pipeline": "BTCRFlow",
        "files": files,
    }

    with open(args.manifest, "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    with open(args.md5, "w") as f:
        for entry in files:
            full_path = os.path.join(workdir, entry["path"])
            if os.path.isfile(full_path):
                checksum = md5_file(full_path)
                f.write(f"{checksum}  {entry['path']}\n")

    print(f"Delivery manifest: {args.manifest}")
    print(f"MD5 checksums: {args.md5}")
    print(f"Files tracked: {len(files)}")


if __name__ == "__main__":
    main()
