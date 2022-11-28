rule NRL_salm:
    input:
        salm_qc = "sample_sheets/QCpass/salmonella.tsv",
        ssheet = "sample_sheets/aquamis_samples.tsv",
    output:
        tab = "exports/NRL_salm/md5_cvua-rrw.txt",
        fastq = directory("exports/NRL_salm"),
        tmp = temp("exports/NRL_salm/filter_table.txt"),
        patterns = temp("exports/NRL_salm/pass_samples.txt"),
    message:
        "Creating exports for NRL Salm"
    log:
        "logs/NRL_salm.log"
    shell:
        """
        exec 2> {log}
        # Get samples id from autoqc and filter aquamis sample sheet
        # then copy to folder
        cut -f 1 {input.salm_qc} | tail -n +2 | tr ' ' \\t > {output.patterns}
        grep -F -w -f {output.patterns} {input.ssheet} > {output.tmp} 
        mkdir -p {output.fastq}
        cp -t $(realpath {output.fastq}) $(awk '{{print $2, $3}}' {output.tmp})
        cd {output.fastq}
        md5sum *.fastq.gz > md5_cvua-rrw.txt
        """
        