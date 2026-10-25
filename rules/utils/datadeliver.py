"""
BTCRFlow target generation (DataDeliver).
Dynamically builds rule all targets based on module switches.
"""
from common import module_enabled, parameter_enabled


MODULE_DEPENDENCIES = {
    "mixcr": ["qc_clean", "trimming"],
    "vdj_analysis": ["mixcr"],
}


def DataDeliver(config, SAMPLES, BCR_SAMPLES, TCR_SAMPLES):
    """Generate target file list for rule all based on enabled modules."""
    targets = []
    all_samples = SAMPLES.index.tolist()

    # Always: seq_preprocessor + MD5 verification (hard gate)
    targets.extend(_convert_md5_targets(all_samples))

    # QC module
    if module_enabled(config, "qc_clean"):
        targets.extend(_qc_targets(all_samples))

    # Trimming module
    if module_enabled(config, "trimming"):
        targets.extend(_trimming_targets(all_samples))

    # MiXCR module
    if module_enabled(config, "mixcr"):
        targets.extend(_mixcr_targets(all_samples, config))

    # VDJ analysis module
    if module_enabled(config, "vdj_analysis"):
        targets.extend(_vdj_targets(all_samples, BCR_SAMPLES, TCR_SAMPLES, config))

    # Delivery
    if module_enabled(config, "deliver"):
        targets.extend(_deliver_targets(config))

    # Report
    if module_enabled(config, "report"):
        targets.extend(_report_targets(config))

    return targets


def _convert_md5_targets(samples):
    """seq_preprocessor + json_md5_verifier outputs (always required)."""
    targets = [
        "00.raw_data/link_dir/raw_data_md5.json",
        "01.qc/md5_check.tsv",
    ]
    for s in samples:
        targets.append(f"00.raw_data/link_dir/{s}/{s}_R1.fq.gz")
        targets.append(f"00.raw_data/link_dir/{s}/{s}_R2.fq.gz")
    return targets


def _qc_targets(samples):
    targets = ["01.qc/multiqc/multiqc.html"]
    return targets


def _trimming_targets(samples):
    targets = []
    for s in samples:
        targets.append(f"02.trimmed/{s}_trimmed_R1.fastq.gz")
        targets.append(f"02.trimmed/{s}_trimmed_R2.fastq.gz")
    return targets


def _mixcr_targets(samples, config):
    targets = []
    for s in samples:
        targets.append(f"03.mixcr/{s}.clns")
    if config.get("mixcr_params", {}).get("postanalysis", True):
        targets.append("03.mixcr/postanalysis.json.gz")
    return targets


def _vdj_targets(samples, bcr_samples, tcr_samples, config):
    targets = []
    vdj_params = config.get("vdjtools_params", {})
    if vdj_params.get("downsample", True):
        for s in tcr_samples:
            targets.append(f"04.vdj_analysis/downsampled/{s}.TRA.downsampled.tsv")
            targets.append(f"04.vdj_analysis/downsampled/{s}.TRB.downsampled.tsv")
        for s in bcr_samples:
            targets.append(f"04.vdj_analysis/downsampled/{s}.IGH.downsampled.tsv")
    if vdj_params.get("diversity", True):
        for s in samples:
            targets.append(f"04.vdj_analysis/diversity/{s}.diversity.done")
    return targets


def _deliver_targets(config):
    return [
        "data_deliver/delivery_manifest.json",
        "data_deliver/delivery_manifest.md5",
    ]


def _report_targets(config):
    return ["05.report/project_summary.json"]
