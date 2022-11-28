import os

rule create_sample_sheet:
    output:
        sample_sheet = "sample_sheets/aquamis_samples.tsv",
    params:
        fastq_folder = config['fastq_folder']
    message:
        "Creating AQUAMIS sample sheet"
    conda:
        "../envs/aquamis.yaml"
    log:
        "logs/aquamis_sample_sheet.log"
    shell:
        """
        exec 2> {log}
        
        create_sampleSheet.sh \
            --mode illumina \
            --fastxDir {params.fastq_folder} \
            --outDir {params.fastq_folder}
        
        mv {params.fastq_folder}/samples.tsv {output.sample_sheet}
        """


rule run_aquamis:
    input:
        sample_sheet = "sample_sheets/aquamis_samples.tsv",
    output:
        report = "AQUAMIS/reports/summary_report.tsv",
        workdir = directory("AQUAMIS"),
    message:
        "Running AQUAMIS"
    conda:
        "../envs/aquamis.yaml"
    threads:
        config['threads'],
    params:
        name = os.path.basename(config['fastq_folder']),
        conda_prefix = config['conda_prefix'],
        min_trim_length = config['min_trim_length'],
        mash_kmer = config['mash_kmer'],
        mash_sketch = config['mash_sketch'],
        bracken_read_length = config['bracken_read_length'],
        shovill_ram = config['shovill_ram'],
        shovill_depth = config['shovill_depth'],
    log:
        "logs/aquamis_run.log"
    shell:
        """
        exec 2> {log}
        
        aquamis \
            --sample_list {input.sample_sheet} \
            --working_directory {output.workdir} \
            --run_name {params.name} \
            --min_trimmed_length {params.min_trim_length} \
            --mash_kmersize {params.mash_kmer} \
            --mash_sketchsize {params.mash_sketch} \
            --read_length {params.bracken_read_length} \
            --shovill_ram {params.shovill_ram} \
            --shovill_depth {params.shovill_depth} \
            --assembler "spades" \
            --threads {threads} \
            --remove_temp \
            --fix_fails 
        """