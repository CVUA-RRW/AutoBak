# First copy assemblies to db folders
# Generate new sample sheet
# Then calculate profiles 

rule copy_lis:
    input:
        ssheet = "sample_sheets/QCpass/listeria.tsv",
    output:
        flag = touch("flags/lis_copyflag"),
    params:
        assemblies = config['lis_assemblies'],
    log:
        "logs/copy_lis_assemblies.log"
    conda:
        "../envs/aquamis.yaml"
    shell:
        """
        exec 2> {log}
        
        cp -t $(realpath {params.assemblies}) $(cut -f 2 {input.ssheet} | tail -n +2)
        
        create_sampleSheet.sh \
            --mode assembly \
            --fastxDir $(realpath {params.assemblies}) \
            --outDir $(realpath {params.assemblies}) \
            --force
        """


rule copy_salm:
    input:
        ssheet = "sample_sheets/QCpass/salmonella.tsv",
    output:
        flag = touch("flags/salm_copyflag"),
    params:
        assemblies = config['salm_assemblies'],
    log:
        "logs/copy_salm_assemblies.log"
    conda:
        "../envs/aquamis.yaml"
    shell:
        """
        exec 2> {log}
        
        cp -t $(realpath {params.assemblies}) $(cut -f 2 {input.ssheet} | tail -n +2)
        
        create_sampleSheet.sh \
            --mode assembly \
            --fastxDir $(realpath {params.assemblies}) \
            --outDir $(realpath {params.assemblies}) \
            --force
        """


rule copy_campy:
    input:
        ssheet = "sample_sheets/QCpass/campylobacter.tsv",
    output:
        flag = touch("flags/campy_copyflag"),
    params:
        assemblies = config['campy_assemblies'],
    log:
        "logs/copy_campy_assemblies.log"
    conda:
        "../envs/aquamis.yaml"
    shell:
        """
        exec 2> {log}
        
        cp -t $(realpath {params.assemblies}) $(cut -f 2 {input.ssheet} | tail -n +2)
        
        create_sampleSheet.sh \
            --mode assembly \
            --fastxDir $(realpath {params.assemblies}) \
            --outDir $(realpath {params.assemblies}) \
            --force
        """


rule copy_coli:
    input:
        ssheet = "sample_sheets/QCpass/escherichia.tsv",
    output:
        flag = touch("flags/coli_copyflag"),
    params:
        assemblies = config['coli_assemblies'],
    log:
        "logs/copy_coli_assemblies.log"
    conda:
        "../envs/aquamis.yaml"
    shell:
        """
        exec 2> {log}
        
        cp -t $(realpath {params.assemblies}) $(cut -f 2 {input.ssheet} | tail -n +2)
        
        create_sampleSheet.sh \
            --mode assembly \
            --fastxDir $(realpath {params.assemblies}) \
            --outDir $(realpath {params.assemblies}) \
            --force
        """


rule merge_lis:
    input:
        flag = "flags/lis_copyflag",
    output:
        flag = touch("flags/Listeria_flag"), 
    params:
        ssheet = f"{config['lis_assemblies']}/samples.tsv", # Species spec
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
        merge_db = os.path.split(os.path.split(config['lis_cgMLST'])[0])[0], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_listeria.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {params.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {params.ssheet} \
            --working_directory {params.merge_db} \
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
            {params.remove_frameshifts} 
        fi
        """


rule merge_salm:
    input:
        flag = "flags/salm_copyflag",
    output:
        flag = touch("flags/Salmonella_flag"), 
    params:
        ssheet = f"{config['salm_assemblies']}/samples.tsv", # Species spec
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
        merge_db = os.path.split(os.path.split(config['salm_cgMLST'])[0])[0], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_salmonella.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {params.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {params.ssheet} \
            --working_directory {params.merge_db} \
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
            {params.remove_frameshifts} 
        fi
        """


rule merge_campy:
    input:
        flag = "flags/campy_copyflag",
    output:
        flag = touch("flags/Campylobacter_flag"), 
    params:
        ssheet = f"{config['campy_assemblies']}/samples.tsv", # Species spec
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
        merge_db = os.path.split(os.path.split(config['campy_cgMLST'])[0])[0], # Species spec
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_campy.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {params.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {params.ssheet} \
            --working_directory {params.merge_db} \
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
            {params.remove_frameshifts} 
        fi
        """


rule merge_coli:
    input:
        flag = "flags/coli_copyflag",
    output:
        flag = touch("flags/Escherichia_flag"), 
    params:
        ssheet = f"{config['coli_assemblies']}/samples.tsv", # Species spec
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
        merge_db = os.path.split(os.path.split(config['coli_cgMLST'])[0])[0], # Species spec    
    threads:
        config['threads'],
    log:
        "logs/chewie_merge_coli.log"
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        """
        exec 2> {log}
        
        if [ $(wc -l < {params.ssheet}) -gt 1 ] ; then
        chewiesnake \
            --threads {threads} \
            --sample_list {params.ssheet} \
            --working_directory {params.merge_db} \
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
            {params.remove_frameshifts} 
        fi
        """


rule flag:
    input:
        "flags/Escherichia_flag",
        "flags/Campylobacter_flag",
        "flags/Listeria_flag",
        "flags/Salmonella_flag",
    output:
        flag = "flags/merge_flag",
    shell:
        "touch {output.flag}"
