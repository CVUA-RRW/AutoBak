
rule precluster:
    input:
        ssheet = "samples_sheets/QCpass_assemblies.tsv",
    output:
        report = "ChewieSnake/report"
    params:
        # sch_listeria = 
        # sch_salmonella = 
        # sch_campylobacter = 
        # sch_escherichia = 
        # pro_listeria =
        # pro_salmonella = 
        # pro_campylobacter = 
        # pro_escherichia = 
    conda:
        "../envs/chewiesnake.yaml"
    shell:
        "touch {output.report}"