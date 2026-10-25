# 08.report.smk
# Lightweight project summary report.


rule project_summary:
    """Generate project_summary.json with run metadata."""
    input:
        mixcr=expand("03.mixcr/{sample}.clns", sample=ALL_SAMPLES),
        multiqc="01.qc/multiqc/multiqc.html",
    output:
        json="05.report/project_summary.json",
    params:
        project_name=config.get("project_name", "btcr_project"),
        species=config.get("species", "Homo_sapiens"),
        genome_version=config.get("Genome_Version", "hg38"),
        n_bcr=len(BCR_SAMPLES),
        n_tcr=len(TCR_SAMPLES),
    log:
        "logs/08.report/project_summary.log",
    threads: 1
    run:
        import json
        import os

        summary = {
            "project_name": params.project_name,
            "species": params.species,
            "genome_version": params.genome_version,
            "total_samples": len(ALL_SAMPLES),
            "bcr_samples": params.n_bcr,
            "tcr_samples": params.n_tcr,
            "sample_list": ALL_SAMPLES,
            "bcr_sample_list": BCR_SAMPLES,
            "tcr_sample_list": TCR_SAMPLES,
            "pipeline": "BTCRFlow",
            "modules": {
                "qc_clean": config.get("qc_clean", False),
                "trimming": config.get("trimming", False),
                "mixcr": config.get("mixcr", False),
                "vdj_analysis": config.get("vdj_analysis", False),
                "deliver": config.get("deliver", False),
            },
        }

        os.makedirs(os.path.dirname(output.json), exist_ok=True)
        with open(output.json, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
