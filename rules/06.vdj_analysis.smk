# 06.vdj_analysis.smk
# VDJtools downstream: downsample, convert, diversity stats.


def get_chains_for_sample(wildcards):
    """Return chain list based on platform."""
    row = SAMPLES.loc[wildcards.sample]
    platform = get_platform(row)
    if platform == "bcr":
        return ["IGH"]
    elif platform == "tcr":
        return ["TRA", "TRB"]
    return ["TRA", "TRB"]


rule mixcr_downsample:
    """Downsample clones to fixed read count per chain."""
    input:
        clns="03.mixcr/{sample}.clns",
    output:
        tsv="04.vdj_analysis/downsampled/{sample}.{chain}.downsampled.tsv",
    params:
        mixcr=get_software(config, "mixcr", "mixcr"),
        count=config.get("vdjtools_params", {}).get("downsample_count", 1000000),
        outdir="04.vdj_analysis/downsampled",
        prefix=lambda wc: f"04.vdj_analysis/downsampled/{wc.sample}.{wc.chain}",
    log:
        "logs/06.vdj/downsample_{sample}_{chain}.log",
    threads: 4
    conda:
        "../envs/mixcr.yaml"
    shell:
        """
        {params.mixcr} exportClones \
            --chains {wildcards.chain} \
            {input.clns} {output.tsv} 2> {log}
        """


rule vdjtools_convert:
    """Convert MiXCR clone TSV to VDJtools format."""
    input:
        tsv="04.vdj_analysis/downsampled/{sample}.{chain}.downsampled.tsv",
    output:
        done="04.vdj_analysis/converted/{sample}.{chain}.convert.done",
    params:
        java=get_software(config, "java", "java"),
        vdjtools=get_software(config, "vdjtools", "vdjtools-1.2.1.jar"),
        mem=config.get("vdjtools_params", {}).get("java_mem", "40g"),
        outdir=lambda wc: f"04.vdj_analysis/converted/{wc.sample}",
    log:
        "logs/06.vdj/convert_{sample}_{chain}.log",
    threads: 2
    run:
        import os
        os.makedirs(params.outdir, exist_ok=True)
        metadata = os.path.join(params.outdir, "metadata.txt")
        with open(metadata, "w") as f:
            f.write("file_name\tsample_id\n")
            f.write(f"{os.path.abspath(input.tsv)}\t{wildcards.sample}\n")
        shell(
            "{params.java} -Xmx{params.mem} -jar {params.vdjtools} "
            "Convert -S MiXcr -c -m {metadata} {params.outdir} 2> {log}"
        )
        shell("touch {output.done}")


rule vdjtools_diversity:
    """Calculate diversity statistics with VDJtools."""
    input:
        done=expand(
            "04.vdj_analysis/converted/{{sample}}.{chain}.convert.done",
            chain=["TRA", "TRB"],
        ),
    output:
        done="04.vdj_analysis/diversity/{sample}.diversity.done",
    params:
        java=get_software(config, "java", "java"),
        vdjtools=get_software(config, "vdjtools", "vdjtools-1.2.1.jar"),
        mem=config.get("vdjtools_params", {}).get("java_mem", "40g"),
        trials=config.get("vdjtools_params", {}).get("diversity_resample_trials", 1),
        outdir=lambda wc: f"04.vdj_analysis/diversity/{wc.sample}",
    log:
        "logs/06.vdj/diversity_{sample}.log",
    threads: 2
    run:
        import os
        import glob as _glob

        os.makedirs(params.outdir, exist_ok=True)
        sample_dir = f"04.vdj_analysis/converted/{wildcards.sample}"
        converted_files = _glob.glob(os.path.join(sample_dir, "*.txt.gz"))
        if not converted_files:
            converted_files = _glob.glob(os.path.join(sample_dir, "*.txt"))

        metadata = os.path.join(params.outdir, "metadata.txt")
        with open(metadata, "w") as f:
            f.write("file_name\tsample_id\n")
            for cf in converted_files:
                f.write(f"{os.path.abspath(cf)}\t{wildcards.sample}\n")

        output_file = os.path.join(params.outdir, "diversity.txt")
        shell(
            "{params.java} -Xmx{params.mem} -jar {params.vdjtools} "
            "CalcDiversityStats --resample-trials {params.trials} "
            "-m {metadata} {output_file} 2> {log}"
        )
        shell("touch {output.done}")
