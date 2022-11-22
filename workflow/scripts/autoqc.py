import sys


sys.stderr = open(snakemake.log[0], "w")


import json
import pandas as pd


def check_criteria(criteria, sample, translator):
    """
    criteria: a dict describing the criteria
    sample: pandas Series
    translator: dict {json_name: aquamis_name}
    
    returns: status in "pass", "warn" or "fail"
    """
    # Get criteria value for sample
    value = sample[translator[criteria['name']]]
    # Evaluate
    if criteria['type'] == 'range':
        flag = value >= criteria['value'][0] and value <= criteria['value'][1]
    elif criteria['type'] == 'threshold':
        flag = value >= criteria['value']
    elif criteria['type'] == 'categorical':
        flag = value ==criteria['value']
    elif criteria['type'] == 'maximum':
        flag = value <= criteria['value']
    else:
        raise ValueError(f"Unknown Argument type for {criteria['name']}: {criteria['type']}")
    # Categorize
    if flag: 
        return "pass"
    else:
        if criteria['category'] == 'obligatory':
            return "fail"
        else:
            return "warn"


def main(metadata, aquamis, thresholds, report):
    meta = pd.read_csv(metadata, sep="\t")
    res = pd.read_csv(aquamis, sep="\t")
    # Creating Genus entry in results table
    res['Genus'] = res.apply(lambda row: row.loc['Species'].split()[0], axis=1)
    with open(thresholds, "r") as handle:
        thr = json.load(handle)
    series = []
    # Iterate over results
    for sample in res.iterrows():
        # retrieve expected species
        sample_id = sample[1]['Sample_Name']
        expect = meta.loc[meta['sample'] == sample_id, 'species'].values[0]
        # Check all QC criteria
        qc_params = thr[expect]
        checks = {entry['name']: check_criteria(entry, sample[1], thr['aquamis_names']) for entry in qc_params}
        # QC_vote
        if 'fail' in checks.values():
            vote = 'Fail'
        elif 'warn' in checks.values():
            vote = 'Warning'
        else:
            vote = "Pass"
        # Create Series
        serie = {"Sample" : sample_id}
        serie.update({'QC vote': vote})
        serie.update(checks)
        series.append(pd.Series(serie))
    df = pd.concat(series, axis=1, ignore_index=True).T
    df.to_csv(report, sep="\t", index=False)


if __name__ == '__main__':
    main(
        metadata=snakemake.input['metadata'],
        aquamis=snakemake.input['aquamis'],
        thresholds=snakemake.params['thresholds'],
        report=snakemake.output['report']
    )