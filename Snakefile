configfile: "config.yaml"

import os
import pandas as pd

workDir = str(os.getcwd())

# Get a list of sample accessions and a dict of sample: run mapping
try:
    accessions = pd.read_csv(workDir + "/data/accessions.tsv", header = 0, sep = '\t')
    SAMPLES = list(accessions.loc[:,'sample_accession'].unique())
    SAMPLES_RUN_DICT = accessions.groupby('sample_accession')['run_accession'].apply(list).to_dict()
except Exception as e:
    print("Error reading accession list: " + workDir + '/data/accessions.tsv')
    print("The file either does not exist, or does not have sample_accession column, which is required.")
    print(e)
    exit(1)


def getPartition(wildcards, resources):
    # Determine partition for each rule based on resources requested
    for key in resources.keys():
        if 'bmm' in key:
            return 'bmm'
        elif 'med' in key:
            return 'med2'
    if int(resources['mem_mb']) / int(resources['cpus']) > 4000:
        return 'bml'
    else:
        return 'low2'
        
rule all:
    input: expand(workDir + "/Results/temp_fastq/{sample}_{read}.fq", sample = SAMPLES, read = ['R1', 'R2'])

include: "rules/prepareFiles.smk"