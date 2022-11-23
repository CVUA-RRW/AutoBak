import os

rule precluster:
    input:
        ssheet = "samples_sheets/QCpass/{species}.tsv",
    output:
        outdir = directory(infer_cgmlst_dir()),
        report = ,
    params:
        schema = infer_schema(),
        prodigal = infer_prodigal(),
        
        pro_listeria = os.path.join({workflow.basedir}, "prodigal_training_files", "Listeria_monocytogenes.trn"),
        pro_salmonella = os.path.join({workflow.basedir}, "prodigal_training_files", "Salmonella_enterica.trn")
        pro_campylobacter = os.path.join({workflow.basedir}, "prodigal_training_files", "Campylobacter_jejuni.trn")
        pro_escherichia = os.path.join({workflow.basedir}, "prodigal_training_files", "Escherichia_coli.trn")
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        chewiesnake \
            --sample_list {input.ssheet} \
            --working_directory {output.outdir} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold
            --size_threshold
            --max_fraction_missing_loci
            --distance_method
            --clustering_method
            --distance_threshold
            --address_range
            --remove_frameshifts
            --allele_length_threshold
            --frameshift_mode
            # Try to integrate comparison here by checking the config file
            {params.comparison} \
            {params.comparison_db} \
            {params.joining_threshold} \
            
            