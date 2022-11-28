
rule merge_lis:
    input:
        ssheet = "sample_sheets/QCpass/listeria.tsv", # Species spec
    output:
        flag = touch("Listeria_flag"), 
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
        merge_db = config['lis_cgMLST'], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_listeria.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {params.merge_db} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            # --max_fraction_missing_loci {params.max_fraction_missing_loci} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            {params.remove_frameshifts} 
        fi
        """


rule merge_salm:
    input:
        ssheet = "sample_sheets/QCpass/salmonella.tsv", # Species spec
    output:
        flag = touch("Salmonella_flag"), 
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
        merge_db = config['salm_cgMLST'], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_salmonella.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {params.merge_db} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            # --max_fraction_missing_loci {params.max_fraction_missing_loci} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            {params.remove_frameshifts} 
        fi
        """


rule merge_campy:
    input:
        ssheet = "sample_sheets/QCpass/campylobacter.tsv", # Species spec
    output:
        flag = touch("Campylobacter_flag"), 
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
        merge_db = config['campy_cgMLST'], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_campy.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {params.merge_db} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            # --max_fraction_missing_loci {params.max_fraction_missing_loci} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            {params.remove_frameshifts} 
        fi
        """


rule merge_coli:
    input:
        ssheet = "sample_sheets/QCpass/escherichia.tsv", # Species spec
    output:
        flag = touch("Escherichia_flag"), 
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
        merge_db = config['coli_cgMLST'], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_coli.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {input.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {input.ssheet} \
            --working_directory {params.merge_db} \
            --scheme {params.schema} \
            --prodigal {params.prodigal} \
            --bsr_threshold {params.bsr_threshold} \
            --size_threshold {params.size_threshold} \
            # --max_fraction_missing_loci {params.max_fraction_missing_loci} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --address_range {params.address_range} \
            --allele_length_threshold {params.allele_length_threshold} \
            --frameshift_mode {params.frameshift_mode} \
            {params.remove_frameshifts} 
        fi
        """

rule flag:
    input:
        "Escherichia_flag",
        "Campylobacter_flag",
        "Listeria_flag",
        "Salmonella_flag",
    output:
        flag = "merge_flag",
    shell:
        "touch {output.flag}"