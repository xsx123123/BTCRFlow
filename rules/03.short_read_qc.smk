# 03.short_read_qc.smk
# FastQC on standardized reads + MultiQC aggregation.
# Depends on md5_check.tsv as hard gate (data integrity verified before QC).


rule fastqc:
    """Run FastQC on standardized linked FASTQ files (after MD5 verification)."""
    input:
        r1="00.raw_data/link_dir/{sample}/{sample}_R1.fq.gz",
        r2="00.raw_data/link_dir/{sample}/{sample}_R2.fq.gz",
        md5_check="01.qc/md5_check.tsv",
    output:
        r1_zip="01.qc/fastqc/{sample}_R1_fastqc.zip",
        r2_zip="01.qc/fastqc/{sample}_R2_fastqc.zip",
    params:
        outdir="01.qc/fastqc",
        fastqc=get_software(config, "fastqc", "fastqc"),
        threads=config.get("qc_params", {}).get("fastqc_threads", 4),
    log:
        "logs/03.qc/fastqc_{sample}.log",
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        {params.fastqc} --outdir {params.outdir} --threads {params.threads} \
            --nogroup --quiet \
            {input.r1} {input.r2} 2> {log}
        """


rule multiqc:
    """Aggregate FastQC results with MultiQC."""
    input:
        zips=expand(
            "01.qc/fastqc/{sample}_R1_fastqc.zip",
            sample=ALL_SAMPLES,
        ) + expand(
            "01.qc/fastqc/{sample}_R2_fastqc.zip",
            sample=ALL_SAMPLES,
        ),
        md5_check="01.qc/md5_check.tsv",
    output:
        html="01.qc/multiqc/multiqc.html",
    params:
        indir="01.qc/fastqc",
        outdir="01.qc/multiqc",
        multiqc=get_software(config, "multiqc", "multiqc"),
    log:
        "logs/03.qc/multiqc.log",
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        {params.multiqc} {params.indir} \
            --outdir {params.outdir} \
            --force \
            --quiet 2> {log}
        """
