"""
BTCRFlow - Bulk BCR/TCR Repertoire Analysis Pipeline
Entry point following the Flow framework conventions.
"""
import os
import sys
import pandas as pd
from pathlib import Path
from snakemake.utils import validate, min_version

min_version("9.9.0")

# --- Configuration loading (Flow family priority) ---
# 1. Default four-piece config
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
configfile: "config/run_parameter.yaml"
configfile: "config/cluster_config.yaml"

# 2. User analysisyaml (highest priority override)
sys.path.insert(0, os.path.join(os.path.dirname(workflow.basedir), "rules", "utils"))
from common import load_user_config, resolve_reference_paths, validate_genome_version, validate_species

load_user_config(config, cmd_arg_name="analysisyaml")
resolve_reference_paths(config)
validate(config, "schema/config.schema.yaml")
validate_genome_version(config)
validate_species(config)

# --- Working directory ---
workdir: config["workflow"]

# --- Sample table ---
from common import load_samples, load_contrasts

SAMPLES = load_samples(config, required_columns=["sample", "sample_name", "group", "platform", "path"])
CONTRASTS = load_contrasts(config)

BCR_SAMPLES = SAMPLES[SAMPLES["platform"].str.lower().isin(["bcr", "bulk_bcr"])].index.tolist()
TCR_SAMPLES = SAMPLES[SAMPLES["platform"].str.lower().isin(["tcr", "bulk_tcr"])].index.tolist()
ALL_SAMPLES = SAMPLES.index.tolist()

# --- Target generation ---
from datadeliver import DataDeliver

rule all:
    input:
        DataDeliver(config, SAMPLES, BCR_SAMPLES, TCR_SAMPLES)

# --- Rules ---
include: "rules/01.common.smk"
include: "rules/02.file_convert_md5.smk"
include: "rules/03.short_read_qc.smk"
include: "rules/04.trimming.smk"
include: "rules/05.mixcr.smk"
include: "rules/06.vdj_analysis.smk"
include: "rules/07.deliver.smk"
include: "rules/08.report.smk"
