# VariantFrequency
Search for variant in ENA/SRA to determine allele frequency

This is an automated Snakemake workflow to search variant allele frequency in a given genome.

At a minimum, it requires [Snakemake](https://snakemake.readthedocs.io/en/stable/), [Python3](https://www.python.org/download/releases/3.0/) and [Pandas](https://pandas.pydata.org/pandas-docs/stable/index.html). It is also recommended to have [Conda](https://docs.conda.io/en/latest/) installed. Otherwise, you need all software listed [here](https://github.com/SichongP/VariantFrequency/blob/main/software_version.txt)

## Get Started

### Install Conda

If you don't already have Conda, it is highly recommended that you install it. This workflow uses Conda to manage rule environments. Although you can still run it without Conda (as long as you make sure all software listed [here](https://github.com/SichongP/VariantFrequency/blob/main/software_version.txt) and their dependencies are installed).

To install Conda:
```
echo source ~/.bashrc >> ~/.bash_profile
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Log out and log back in again to activate the base conda environment.


### Adjust configurations

The main snakemake workflow is in [`Snakefile`](https://github.com/SichongP/VariantFrequency/blob/main/Snakefile) inside the main directory. Each submodule is stored in [`rules/`](https://github.com/SichongP/VariantFrequency/tree/main/rules).

The workflow starts by downloading reference genome FASTA file. The genome can be specified in [`config.yaml`](https://github.com/SichongP/VariantFrequency/blob/main/config.yaml):
```
genome: hg38
```
Make sure that the genome name is exactly same as in corresponding database (eg. equCab3, hg19, mm10, etc)

Alternatively, if you have a reference genome FASTA file already, you can also put it under `data/` directory with the name `genome.fa`. This will bypass the step to download reference genome.

A list of variants to search is stored in [`data/variants.bed`](https://github.com/SichongP/VariantFrequency/blob/main/data/variants.bed). Modify this list to include the variant positions you want to search.
Keep in mind that `bed` files are 0-based (each chromosome positions starts at 0, the first position on chr1 should be chr1:0-1). This file should have three columns, tab-delimited, with no header.

A list of accessions is stored in [`data/accessions.tsv`](https://github.com/SichongP/VariantFrequency/blob/main/data/accessions.tsv). This is a tab-delimited text file with at least the following columns:
```
- sample_accession
- run_accession
```
`run_accession` allows us to download fastq files for each run while `sample_accession` allows us to combine different runs of same sample together. [`ENA_curl.sh`](https://github.com/SichongP/VariantFrequency/blob/main/ENA_curl.sh) contains an example ENA query you can use. Modify the query string to match your criteria.
To generate an accession file, simply run: `./ENA_curl.sh > data/accessions.tsv`

### Running the workflow

To start the Snakemake workflow on a local machine or a dedicated server with Conda:
```
snakemake --cores n_cores --use-conda
```
Without Conda:
```
snakemake --cores n_cores
```
`n_cores` is the number of threads Snakemake is allowed to use for this workflow.

After the workflow successfully completes, an output file `output.vcf` can be found in `Results/VCF/` #Update this after adding rule to calculate AF

### Running on HPC

You can also run this workflow on an HPC cluster. Read [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for more details.
Most rules in this workflow have some common resources defined (memory, cpus, and time) under `resources` directive. You can access them using these keywords:
```
mem_mb: memory in MB
cpus: threads
time_min: time in minutes
```
Additionally, a `partition` keyword is defined in `params` directive. Since this is usually cluster-specific, you will likely want to modify the [`getPartition`](https://github.com/SichongP/VariantFrequency/blob/41a3b157c93cf5015f3ade5648680469c3f9ee8e/Snakefile#L20-L30) function defined in the main `Snakefile`.
The current configuration should work out-of-box on UC Davis Farm cluster.

### Performance

1. It's taking too long  
This workflow downloads raw fastq files available on SRA/ENA database and aligns them to a list of reference sequences containing 1kb up- and down-stream of each variant site. The main bottle neck is at the download step.
If storage space is not a concern, it is advisable to keep these fastq files so future runs can skip the download step. (The default behaviour is to remove them after alignment to save space.)
To do this:
delete the `temp()` directive surrounding the below lines in [`rules/prepareFiles.smk`](https://github.com/SichongP/VariantFrequency/blob/main/rules/prepareFiles.smk) (line 90-91)
```
        r1 = temp(workDir + "/Results/temp_fastq/{sample}_R1.fq"),
        r2 = temp(workDir + "/Results/temp_fastq/{sample}_R2.fq")
```

2. The Snakemake hangs for a long time with no output
If the accession list is too long, Snakemake can take a long time to build DAG for the workflow. This is mostly due to ecessive queries on the drive for output files. You can improve
performance by running the workflow in batches:

```
snakemake --cores n_cores --use-conda --batch genotypeGVCF=1/3
```
will separate the workflow into 3 batches and run the first batch. Once that is completed, you can then run:
```
snakemake --cores n_cores --use-conda --batch genotypeGVCF=2/3
```
and finally:
```
snakemake --cores n_cores --use-conda --batch genotypeGVCF=3/3
```


