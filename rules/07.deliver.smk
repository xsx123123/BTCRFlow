# 07.deliver.smk
# Delivery manifest generation and file delivery.


rule generate_delivery_manifest:
    """Generate delivery manifest JSON and MD5 checksums."""
    input:
        mixcr=expand("03.mixcr/{sample}.clns", sample=ALL_SAMPLES),
        multiqc="01.qc/multiqc/multiqc.html",
    output:
        manifest="data_deliver/delivery_manifest.json",
        md5="data_deliver/delivery_manifest.md5",
    params:
        deliver_script=get_software(config, "btcrflow_deliver", "../scripts/btcrflow_deliver.py"),
        data_deliver=config.get("data_deliver", "data_deliver"),
        project_name=config.get("project_name", "btcr_project"),
    log:
        "logs/07.deliver/delivery.log",
    threads: 1
    shell:
        """
        python {params.deliver_script} \
            --project-name {params.project_name} \
            --output-dir {params.data_deliver} \
            --manifest {output.manifest} \
            --md5 {output.md5} \
            2> {log}
        """
