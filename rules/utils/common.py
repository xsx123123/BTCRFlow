"""
BTCRFlow common utilities.
Configuration loading, sample table parsing, and validation helpers.
"""
import os
import sys
import yaml
import pandas as pd
from pathlib import Path


def load_user_config(config, cmd_arg_name="analysisyaml"):
    """Load user analysisyaml and overlay onto default config (highest priority)."""
    analysisyaml = config.get(cmd_arg_name)
    if not analysisyaml:
        return
    if not os.path.isfile(analysisyaml):
        raise FileNotFoundError(f"analysisyaml not found: {analysisyaml}")
    with open(analysisyaml) as f:
        user_cfg = yaml.safe_load(f) or {}
    _deep_update(config, user_cfg)


def _deep_update(base, override):
    """Recursively merge override into base dict."""
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            _deep_update(base[key], value)
        else:
            base[key] = value


def resolve_reference_paths(config, index_keys=None):
    """Resolve reference_path into config if reference block exists."""
    ref_path = config.get("reference_path")
    if not ref_path:
        return
    if index_keys is None:
        index_keys = ["mixcr_species", "vdjtools_jar"]
    ref_block = config.get("reference", {})
    for key in index_keys:
        if key in ref_block:
            config.setdefault("resolved_references", {})[key] = ref_block[key]


def validate_genome_version(config):
    """Validate genome version is set."""
    gv = config.get("Genome_Version")
    if not gv:
        raise ValueError("Genome_Version must be specified in config or analysisyaml")


def validate_species(config):
    """Validate species is supported."""
    supported = ["Homo_sapiens", "Mus_musculus"]
    sp = config.get("species")
    if sp not in supported:
        raise ValueError(f"species '{sp}' not in supported list: {supported}")


def load_samples(config, required_columns=None):
    """Load and validate sample CSV."""
    if required_columns is None:
        required_columns = ["sample", "sample_name", "group", "platform", "path"]
    sample_csv = config.get("sample_csv")
    if not sample_csv:
        raise ValueError("sample_csv not specified")
    if not os.path.isfile(sample_csv):
        raise FileNotFoundError(f"sample_csv not found: {sample_csv}")
    samples = pd.read_csv(sample_csv, comment="#", dtype=str).fillna("")
    missing = [c for c in required_columns if c not in samples.columns]
    if missing:
        raise ValueError(f"samples.csv missing required columns: {missing}")
    samples = samples.set_index("sample", drop=False)
    return samples


def load_contrasts(config):
    """Load contrasts CSV (optional)."""
    paired_csv = config.get("paired_csv")
    if not paired_csv or not os.path.isfile(paired_csv):
        return pd.DataFrame()
    return pd.read_csv(paired_csv, comment="#", dtype=str).fillna("")


def get_platform(sample_row):
    """Normalize platform string to bcr/tcr."""
    p = str(sample_row.get("platform", "")).lower().strip()
    if p in ("bcr", "bulk_bcr", "bulk-bcr"):
        return "bcr"
    if p in ("tcr", "bulk_tcr", "bulk-tcr"):
        return "tcr"
    return p


def get_mixcr_preset(config, platform):
    """Return MiXCR preset name based on platform."""
    params = config.get("mixcr_params", {})
    if platform == "bcr":
        return params.get("bcr_preset", "takara-human-rna-bcr-umi-smarter")
    elif platform == "tcr":
        return params.get("tcr_preset", "generic-amplicon")
    return params.get("tcr_preset", "generic-amplicon")


def get_sample_path(sample_row):
    """Get raw data path for a sample."""
    return sample_row.get("path", "")


def get_sample_index(sample_row):
    """Get file index/prefix for a sample."""
    idx = sample_row.get("index", "")
    if not idx:
        return sample_row.get("sample", "")
    return idx


def module_enabled(config, module_name):
    """Check if a module is enabled."""
    return config.get(module_name, False)


def parameter_enabled(config, param_path):
    """Check nested parameter run flag, e.g. 'trimming.run'."""
    parts = param_path.split(".")
    val = config.get("parameter", {})
    for p in parts:
        if isinstance(val, dict):
            val = val.get(p, {})
        else:
            return False
    return bool(val) if not isinstance(val, dict) else val.get("run", False)
