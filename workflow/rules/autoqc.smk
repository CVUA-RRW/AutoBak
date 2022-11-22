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
            --fastxDir {params.fasta_folder} \
            --outDir sample_sheets \
            --force
        
        mv sample_sheets/samples.tsv {output.assemblies}
        """


rule QC_filter_sample_sheet:
    input:
        autoqc = "AutoQC/autoqc_check.tsv",
        assemblies = "sample_sheets/all_assemblies.tsv",
    output:
        ssheet = "samples_sheets/QCpass_assemblies.tsv",
    params:
        keep_warn = config['keep_warn'],
    conda:
        "../envs/pandas.yaml"
    message:
        "Filtering sample sheet"
    log:
       "logs/autoQC_ssheet.log"
    script:
        "../scripts/qc_sample_sheet.py"
        


