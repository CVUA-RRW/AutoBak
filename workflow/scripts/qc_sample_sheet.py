import sys


sys.stderr = open(snakemake.log[0], "w")


import pandas as pd


def main(autoqc, assemblies, ssheet, keep_warn):
    qcvotes = pd.read_csv(autoqc, sep="\t")
    paths = pd.read_csv(assemblies, sep="\t")
    filter = ['Pass']
    if keep_warn:
        filter.append('Warning')
    keep = qcvotes.loc[qcvotes['QC vote'].isin(filter)]['Sample'].tolist()
    df = paths.loc[paths['sample'].isin(keep)]
    df.to_csv(ssheet, sep="\t", index=False)

if __name__ == '__main__':
    main(
        autoqc=snakemake.input['autoqc'],
        assemblies=snakemake.input['assemblies'],
        keep_warn=snakemake.params['keep_warn'],
        ssheet=snakemake.output['ssheet']
    )