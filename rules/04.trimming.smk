# 04.trimming.smk
# Adapter trimming with cutadapt (Clontech BCR/TCR protocol).
# Input from standardized link_dir paths (post MD5 verification).


rule cutadapt_trim:
    """Trim adapters and hard-trim 5' end for Clontech SMARTer protocol."""
    input:
        r1="00.raw_data/link_dir/{sample}/{sample}_R1.fq.gz",
        r2="00.raw_data/link_dir/{sample}/{sample}_R2.fq.gz",
        md5_check="01.qc/md5_check.tsv",
    output:
        r1="02.trimmed/{sample}_trimmed_R1.fastq.gz",
        r2="02.trimmed/{sample}_trimmed_R2.fastq.gz",
        json="02.trimmed/{sample}_trim_report.json",
    params:
        cutadapt=get_software(config, "cutadapt", "cutadapt"),
        hard_trim=config.get("trimming_params", {}).get("hard_trim_5prime", 30),
        quality=config.get("trimming_params", {}).get("quality", "25,25"),
        error_rate=config.get("trimming_params", {}).get("error_rate", 0.1),
        min_length=config.get("trimming_params", {}).get("min_length", 20),
        fwd_adapter=config.get("trimming_params", {}).get(
            "fwd_adapter", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        ),
        rev_adapter=config.get("trimming_params", {}).get(
            "rev_adapter", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        ),
    log:
        "logs/04.trimming/{sample}.log",
    threads: 8
    conda:
        "../envs/trimming.yaml"
    shell:
        """
        {params.cutadapt} \
            -U {params.hard_trim} \
            --trim-n \
            -q {params.quality} \
            -e {params.error_rate} \
            -m {params.min_length} \
            -a {params.fwd_adapter} \
            -A {params.rev_adapter} \
            -j {threads} \
            --json={output.json} \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} 2> {log}
        """
