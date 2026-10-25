"""
BTCRFlow tools helpers.
Wrapper functions for constructing shell commands.
"""
import os


def get_software(config, key, default=None):
    """Get software path from config."""
    sw = config.get("software", {})
    return sw.get(key, default or key)


def mixcr_analyze_cmd(config, preset, input_r1, input_r2, output_prefix, threads=12):
    """Build MiXCR analyze command."""
    mixcr = get_software(config, "mixcr", "mixcr")
    params = config.get("mixcr_params", {})
    species = params.get("species", "hsa")

    cmd_parts = [mixcr, "analyze", f'"{preset}"']
    cmd_parts.append(f"--species {species}")
    if params.get("rna", True):
        cmd_parts.append("--rna")
    if params.get("rigid_left_alignment_boundary", True):
        cmd_parts.append("--rigid-left-alignment-boundary")
    boundary = params.get("floating_right_alignment_boundary", "C")
    if boundary:
        cmd_parts.append(f"--floating-right-alignment-boundary {boundary}")
    cmd_parts.append(f"-t {threads}")
    cmd_parts.append(f"{input_r1} {input_r2} {output_prefix}")
    return " \\\n    ".join(cmd_parts)


def mixcr_postanalysis_cmd(config, clns_files, output_file):
    """Build MiXCR postanalysis command."""
    mixcr = get_software(config, "mixcr", "mixcr")
    files_str = " ".join(clns_files)
    return (
        f"{mixcr} postanalysis individual "
        f"--default-downsampling none "
        f"--default-weight-function none "
        f"{files_str} {output_file}"
    )


def mixcr_downsample_cmd(config, input_clns, chain, output_prefix, count=1000000):
    """Build MiXCR downsample command."""
    mixcr = get_software(config, "mixcr", "mixcr")
    return (
        f"{mixcr} downsample "
        f"--chains {chain} "
        f"--downsampling count-read-fixed-{count} "
        f"--suffix downsampled "
        f"{input_clns} {output_prefix}"
    )


def mixcr_export_clones_cmd(config, input_clns, chain, output_tsv):
    """Build MiXCR exportClones command."""
    mixcr = get_software(config, "mixcr", "mixcr")
    return (
        f"{mixcr} exportClones "
        f"--chains {chain} "
        f"{input_clns} {output_tsv}"
    )


def vdjtools_convert_cmd(config, input_tsv, output_dir, metadata_file):
    """Build VDJtools Convert command."""
    java = get_software(config, "java", "java")
    vdjtools = get_software(config, "vdjtools", "vdjtools-1.2.1.jar")
    mem = config.get("vdjtools_params", {}).get("java_mem", "40g")
    return (
        f"{java} -Xmx{mem} -jar {vdjtools} Convert "
        f"-S MiXcr -c -m {metadata_file} {output_dir}"
    )


def vdjtools_diversity_cmd(config, metadata_file, output_file):
    """Build VDJtools CalcDiversityStats command."""
    java = get_software(config, "java", "java")
    vdjtools = get_software(config, "vdjtools", "vdjtools-1.2.1.jar")
    mem = config.get("vdjtools_params", {}).get("java_mem", "40g")
    trials = config.get("vdjtools_params", {}).get("diversity_resample_trials", 1)
    return (
        f"{java} -Xmx{mem} -jar {vdjtools} CalcDiversityStats "
        f"--resample-trials {trials} "
        f"-m {metadata_file} {output_file}"
    )


def cutadapt_cmd(config, input_r1, input_r2, output_r1, output_r2, threads=4):
    """Build cutadapt command for Clontech BCR/TCR trimming."""
    cutadapt = get_software(config, "cutadapt", "cutadapt")
    params = config.get("trimming_params", {})
    hard_trim = params.get("hard_trim_5prime", 30)
    quality = params.get("quality", "25,25")
    error_rate = params.get("error_rate", 0.1)
    min_length = params.get("min_length", 20)
    fwd_adapter = params.get("fwd_adapter", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
    rev_adapter = params.get("rev_adapter", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")

    cmd_parts = [cutadapt]
    cmd_parts.append(f"-U {hard_trim}")
    if params.get("trim_n", True):
        cmd_parts.append("--trim-n")
    cmd_parts.append(f"-q {quality}")
    cmd_parts.append(f"-e {error_rate}")
    cmd_parts.append(f"-m {min_length}")
    cmd_parts.append(f"-a {fwd_adapter}")
    cmd_parts.append(f"-A {rev_adapter}")
    cmd_parts.append(f"-j {threads}")
    cmd_parts.append(f"-o {output_r1} -p {output_r2}")
    cmd_parts.append(f"{input_r1} {input_r2}")
    return " \\\n    ".join(cmd_parts)
