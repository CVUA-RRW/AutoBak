import os 


rule quality_control:
    input:
        metadata = "Metadata/metadata.tsv",
        aquamis = "AQUAMIS/reports/summary_report.tsv",
    output:
        report = "AutoQC/autoqc_check.tsv",
    params:
        thresholds = os.path.join(workflow.basedir, "..", "data", "AutoQC_thresholds.json"),
    message:
        "Running AutoQC"
    conda:
        "../envs/pandas.yaml"
    log:
       "logs/autoQC.log"
    script:
        "../scripts/autoqc.py"


rule assembly_sample_sheet:
    input:
        aquamis = "AQUAMIS/reports/summary_report.tsv",
    output:
        assemblies = "sample_sheets/all_assemblies.tsv",
    params:
        fasta_folder = "AQUAMIS/Assembly/assembly"
    message:
        "Creating assembly sample sheet"
    conda:
        "../envs/aquamis.yaml"
    log:
        "logs/assmebly_sample_sheet.log"
    shell:
        """
        exec 2> {log}
        
        create_sampleSheet.sh \
            --mode assembly \
            --fastxDir $(realpath {params.fasta_folder}) \
            --outDir sample_sheets \
            --force
        
        mv sample_sheets/samples.tsv {output.assemblies}
        """


rule QC_filter_sample_sheet:
    input:
        autoqc = "AutoQC/autoqc_check.tsv",
        assemblies = "sample_sheets/all_assemblies.tsv",
        metadata = "Metadata/metadata.tsv",
    output:
        out_lis = "sample_sheets/QCpass/listeria.tsv",
        out_salm = "sample_sheets/QCpass/salmonella.tsv",
        out_campy = "sample_sheets/QCpass/campylobacter.tsv",
        out_ecoli = "sample_sheets/QCpass/escherichia.tsv",
    params:
        keep_warn = config['keep_warn'],
        thresholds = os.path.join(workflow.basedir, "..", "data", "AutoQC_thresholds.json"),
    conda:
        "../envs/pandas.yaml"
    message:
        "Filtering sample sheet"
    log:
       "logs/autoQC_ssheet.log"
    script:
        "../scripts/qc_sample_sheet.py"
        


rule autoqc_report:
    input:
        autoqc = "AutoQC/autoqc_check.tsv",
    output:
        autoqc = "AutoQC/autoqc.html",
    params:
        workdir = config['workdir'],
        version = version,
    conda:
        "../envs/aquamis.yaml"
    message:
        "Writting AutoQC report"
    log:
        "logs/write_report.log"
    script:
        ".../scripts/write_autoqc_report.Rmd"

    