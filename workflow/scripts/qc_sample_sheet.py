import sys


sys.stderr = open(snakemake.log[0], "w")


import os
import pandas as pd


def main(autoqc, assemblies, metadata, keep_warn, outdir):
    qcvotes = pd.read_csv(autoqc, sep="\t")
    paths = pd.read_csv(assemblies, sep="\t")
    metadata = pd.read_csv(metadata, sep='\t')
    filter = ['Pass']
    if keep_warn:
        filter.append('Warning')
    # Filter by QC criteria
    keep = qcvotes.loc[qcvotes['QC vote'].isin(filter)]['Sample'].tolist()
    # Split by species and produce one sample sheet per species
    for species in metadata['species'].unique():
        df = paths.loc[(paths['sample'].isin(keep)) & (paths['sample']==species)]
        df.to_csv(os.path.join(outdir, f"{species}.tsv"), sep="\t", index=False)

if __name__ == '__main__':
    main(
        autoqc=snakemake.input['autoqc'],
        assemblies=snakemake.input['assemblies'],
        metadata=snakemake.input['metadata'],
        keep_warn=snakemake.params['keep_warn'],
        outdir=snakemake.output['outdir']
    )