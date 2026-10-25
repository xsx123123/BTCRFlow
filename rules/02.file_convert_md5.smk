# 02.file_convert_md5.smk
# FASTQ naming standardization (seq_preprocessor) + MD5 verification (json_md5_verifier).
# Follows RNAFlow three-step pattern: standardize -> verify -> gate for downstream.


def get_all_input_dirs(sample_keys, config):
    """Collect unique raw data directories from sample table and config."""
    dirs = set()
    for s in sample_keys:
        row = SAMPLES.loc[s]
        p = get_sample_path(row)
        if p:
            dirs.add(p)
    for p in config.get("raw_data_path", []):
        dirs.add(p)
    return sorted(dirs)


rule seq_preprocessor:
    """Standardize FASTQ naming and generate raw_data_md5.json."""
    input:
        dirs=get_all_input_dirs(ALL_SAMPLES, config),
    output:
        json_report="00.raw_data/link_dir/raw_data_md5.json",
        link_dir=directory("00.raw_data/link_dir"),
    params:
        seq_preprocessor=get_software(config, "seq_preprocessor"),
        md5_name=config.get("raw_data", {}).get("md5", "md5.txt"),
        library_type=config.get("raw_data", {}).get("library_type", "short-read"),
        raw_data_path=lambda wc: " ".join(get_all_input_dirs(ALL_SAMPLES, config)),
    log:
        "logs/02.convert_md5/seq_preprocessor.log",
    threads: 4
    shell:
        """
        {params.seq_preprocessor} \
            -i {params.raw_data_path} \
            -o {output.link_dir} \
            --md5-name {params.md5_name} \
            --library-type {params.library_type} \
            --json-report {output.json_report} \
            2> {log}
        """


rule check_md5:
    """Verify MD5 checksums from raw_data_md5.json (hard gate for downstream)."""
    input:
        json_report="00.raw_data/link_dir/raw_data_md5.json",
    output:
        tsv="01.qc/md5_check.tsv",
    params:
        verifier=get_software(config, "json_md5_verifier"),
        base_dir="00.raw_data/link_dir",
    log:
        "logs/02.convert_md5/check_md5.log",
    threads: 40
    shell:
        """
        {params.verifier} \
            -t {threads} \
            -i {input.json_report} \
            -b {params.base_dir} \
            -o {output.tsv} \
            --log-file {log}
        """
