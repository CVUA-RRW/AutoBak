rule copy_metadata:
    input:
        meta = config['metadata'],
    output: 
        metadata = "Metadata/metadata.tsv",
    message:
        "Copying metadata"
    shell:
        "cp {input.meta} {output.metadata}"
