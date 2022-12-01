
rule compare_lis:
    input:
        ssheet = "sample_sheets/QCpass/listeria.tsv", # Species spec
    output:
        outdir = directory("cgMLST_Listeria"), # Species spec
        report = "cgMLST_Listeria/reports/cgmlst_comparison_report.html", # Species spec
    params:
        schema = config['lis_schema'], # Species spec
        prodigal = config['lis_prodigal'], # Species spec
        bsr_threshold = config['bsr_threshold'],
        size_threshold = config['size_threshold'],
        max_fraction_missing_loci = config['max_fraction_missing_loci'],
        distance_method = config['distance_method'],
        clustering_method = config['clustering_method'],
        distance_threshold = config['lis_dist'], # Species spec
        address_range = config['address_range'],
        remove_frameshifts = "--remove_frameshift" if config['remove_frameshifts'] else "",
        allele_length_threshold = config['allele_length_threshold'],
        frameshift_mode = config['frameshift_mode'],
        compare_db = config['lis_cgMLST'], # Species spec
        joining_threshold = config['joining_threshold'],
    threads:
        config['threads'],
    log:
        "logs/chewie_run_listeria.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {output.outdir} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            --comparison --comparison_db {params.compare_db} --report \
            --joining_threshold {params.joining_threshold} \
            {params.remove_frameshifts} 
        fi
        """


rule compare_salm:
    input:
        ssheet = "sample_sheets/QCpass/salmonella.tsv", # Species spec
    output:
        outdir = directory("cgMLST_Salmonella"), # Species spec
        report = "cgMLST_Salmonella/reports/cgmlst_comparison_report.html", # Species spec
    params:
        schema = config['salm_schema'], # Species spec
        prodigal = config['salm_prodigal'], # Species spec
        bsr_threshold = config['bsr_threshold'],
        size_threshold = config['size_threshold'],
        max_fraction_missing_loci = config['max_fraction_missing_loci'],
        distance_method = config['distance_method'],
        clustering_method = config['clustering_method'],
        distance_threshold = config['salm_dist'], # Species spec
        address_range = config['address_range'],
        remove_frameshifts = "--remove_frameshift" if config['remove_frameshifts'] else "",
        allele_length_threshold = config['allele_length_threshold'],
        frameshift_mode = config['frameshift_mode'],
        compare_db = config['salm_cgMLST'], # Species spec
        joining_threshold = config['joining_threshold'],
    threads:
        config['threads'],
    log:
        "logs/chewie_run_salmonella.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {output.outdir} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            --comparison --comparison_db {params.compare_db} --report \
            --joining_threshold {params.joining_threshold} \
            {params.remove_frameshifts} 
        fi
        """


rule compare_campy:
    input:
        ssheet = "sample_sheets/QCpass/campylobacter.tsv", # Species spec
    output:
        outdir = directory("cgMLST_Campylobacter"), # Species spec
        report = "cgMLST_Campylobacter/reports/cgmlst_comparison_report.html", # Species spec
    params:
        schema = config['campy_schema'], # Species spec
        prodigal = config['campy_prodigal'], # Species spec
        bsr_threshold = config['bsr_threshold'],
        size_threshold = config['size_threshold'],
        max_fraction_missing_loci = config['max_fraction_missing_loci'],
        distance_method = config['distance_method'],
        clustering_method = config['clustering_method'],
        distance_threshold = config['campy_dist'], # Species spec
        address_range = config['address_range'],
        remove_frameshifts = "--remove_frameshift" if config['remove_frameshifts'] else "",
        allele_length_threshold = config['allele_length_threshold'],
        frameshift_mode = config['frameshift_mode'],
        compare_db = config['campy_cgMLST'], # Species spec
        joining_threshold = config['joining_threshold'],
    threads:
        config['threads'],
    log:
        "logs/chewie_run_campy.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {output.outdir} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            --comparison --comparison_db {params.compare_db} --report \
            --joining_threshold {params.joining_threshold} \
            {params.remove_frameshifts} 
        fi
        """


rule compare_coli:
    input:
        ssheet = "sample_sheets/QCpass/escherichia.tsv", # Species spec
    output:
        outdir = directory("cgMLST_Escherichia"), # Species spec
        report = "cgMLST_Escherichia/reports/cgmlst_comparison_report.html", # Species spec
    params:
        schema = config['coli_schema'], # Species spec
        prodigal = config['coli_prodigal'], # Species spec
        bsr_threshold = config['bsr_threshold'],
        size_threshold = config['size_threshold'],
        max_fraction_missing_loci = config['max_fraction_missing_loci'],
        distance_method = config['distance_method'],
        clustering_method = config['clustering_method'],
        distance_threshold = config['coli_dist'], # Species spec
        address_range = config['address_range'],
        remove_frameshifts = "--remove_frameshift" if config['remove_frameshifts'] else "",
        allele_length_threshold = config['allele_length_threshold'],
        frameshift_mode = config['frameshift_mode'],
        compare_db = config['coli_cgMLST'], # Species spec
        joining_threshold = config['joining_threshold'],
    threads:
        config['threads'],
    log:
        "logs/chewie_run_coli.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {output.outdir} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            --comparison --comparison_db {params.compare_db} --report \
            --joining_threshold {params.joining_threshold} \
            {params.remove_frameshifts} 
        fi
        """

