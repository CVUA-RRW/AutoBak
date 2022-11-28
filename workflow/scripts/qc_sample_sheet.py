import sys


sys.stderr = open(snakemake.log[0], "w")


import os
import pandas as pd


def main(autoqc, assemblies, metadata, keep_warn, thresholds, out_lis, out_salm, out_campy, out_ecoli):
    qcvotes = pd.read_csv(autoqc, sep="\t")
    paths = pd.read_csv(assemblies, sep="\t")
    metadata = pd.read_csv(metadata, sep='\t')
    filter = ['Pass']
    if keep_warn:
        filter.append('Warning')
    # Filter by QC criteria
    keep = qcvotes.loc[qcvotes['QC vote'].isin(filter)]
    # Split by species and produce one sample sheet per species
    # Listeria
    keep_spc = keep.loc[keep['species'] == "Listeria monocytogenes"]['Sample'].tolist()
    df = paths.loc[(paths['sample'].isin(keep_spc))]
    df.to_csv(os.path.join(out_lis), sep="\t", index=False)
    # Salmonella
    keep_spc = keep.loc[keep['species'] == "Salmonella enterica"]['Sample'].tolist()
    df = paths.loc[(paths['sample'].isin(keep_spc))]
    df.to_csv(os.path.join(out_salm), sep="\t", index=False)
    # Campys
    keep_spc = keep.loc[keep['species'] == "Campylobacter ssp"]['Sample'].tolist()
    df = paths.loc[(paths['sample'].isin(keep_spc))]
    df.to_csv(os.path.join(out_campy), sep="\t", index=False)
    # Ecoli
    keep_spc = keep.loc[keep['species'] == "Escherichia coli"]['Sample'].tolist()
    df = paths.loc[(paths['sample'].isin(keep_spc))]
    df.to_csv(os.path.join(out_ecoli), sep="\t", index=False)


if __name__ == '__main__':
    main(
        autoqc=snakemake.input['autoqc'],
        assemblies=snakemake.input['assemblies'],
        metadata=snakemake.input['metadata'],
        keep_warn=snakemake.params['keep_warn'],
        thresholds=snakemake.params['thresholds'],
        out_lis=snakemake.output['out_lis'],
        out_salm=snakemake.output['out_salm'],
        out_campy=snakemake.output['out_campy'],
        out_ecoli=snakemake.output['out_ecoli']
    )