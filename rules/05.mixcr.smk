# 05.mixcr.smk
# MiXCR analyze (BCR/TCR preset selection) + postanalysis.
# Input from standardized link_dir paths (post MD5 verification).


def get_mixcr_preset_for_sample(wildcards):
    row = SAMPLES.loc[wildcards.sample]
    platform = get_platform(row)
    return get_mixcr_preset(config, platform)


rule mixcr_analyze:
    """Run MiXCR analyze with platform-specific preset."""
    input:
        r1="00.raw_data/link_dir/{sample}/{sample}_R1.fq.gz",
        r2="00.raw_data/link_dir/{sample}/{sample}_R2.fq.gz",
        md5_check="01.qc/md5_check.tsv",
    output:
        clns="03.mixcr/{sample}.clns",
    params:
        preset=get_mixcr_preset_for_sample,
        mixcr=get_software(config, "mixcr", "mixcr"),
        species=config.get("mixcr_params", {}).get("species", "hsa"),
        outdir="03.mixcr",
        prefix=lambda wc: f"03.mixcr/{wc.sample}",
    log:
        "logs/05.mixcr/{sample}.log",
    threads: lambda: config.get("mixcr_params", {}).get("threads", 12)
    conda:
        "../envs/mixcr.yaml"
    shell:
        """
        {params.mixcr} analyze "{params.preset}" \
            --species {params.species} \
            --rna \
            --rigid-left-alignment-boundary \
            --floating-right-alignment-boundary C \
            -t {threads} \
            {input.r1} {input.r2} {params.prefix} 2> {log}
        """


rule mixcr_postanalysis:
    """Run MiXCR postanalysis on all clone files."""
    input:
        clns=expand("03.mixcr/{sample}.clns", sample=ALL_SAMPLES),
    output:
        json_gz="03.mixcr/postanalysis.json.gz",
    params:
        mixcr=get_software(config, "mixcr", "mixcr"),
    log:
        "logs/05.mixcr/postanalysis.log",
    threads: 4
    conda:
        "../envs/mixcr.yaml"
    shell:
        """
        {params.mixcr} postanalysis individual \
            --default-downsampling none \
            --default-weight-function none \
            {input.clns} {output.json_gz} 2> {log}
        """
