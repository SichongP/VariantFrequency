localrules: getGenomeFASTA, generateRef, generateFASTQLink

def getLink(wildcards, input):
    run, read = wildcards.name.split('_')
    accessions = pd.read_csv(input.Links, header = 0)
    accessions = accessions[accessions['run_accession'] == run]
    if read == '1':
        return accessions['r1_fastq'].iloc[0]
    else:
        return accessions['r2_fastq'].iloc[0]

rule getGenomeFASTA:
    output: fa = workDir + "/data/genome.fa", size = workDir + "/data/chrom.sizes"
    params:
        build = lambda wildcards: config['refBuild'],
        genome = lambda wildcards: config['genome']
    shell:
     """
     if [ "{params.build}" = "UCSC" ]; then
         rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/{params.genome}.fa.gz {output.fa}.gz
         rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/{params.genome}.chrom.sizes {output.size}
         gunzip {output.fa}.gz
     fi
     """
     
rule generateRef:
    input:
        ref = workDir + "/data/genome.fa",
        bed = workDir + "/data/variants.bed",
        chromsize = workDir + "/data/chrom.sizes"
    output:
        fa = workDir + "/data/regionsRef.fa",
        remapVar = workDir + "/data/remappedVariants.bed",
        mapping = workDir + "/data/variantMapping.csv"
    conda: workDir + "/envs/python.yaml"
    script: workDir + "/scripts/generateRef.py"
    
rule generateFASTQLink:
    input: accession = workDir + "/data/accessions.tsv"
    output: outLinks = workDir + "/data/fastqLinks.csv"
    conda: workDir + "/envs/python.yaml"
    script: workDir + "/scripts/getFASTQLinks.py"
    
rule getFASTQ:
    input: Links = workDir + "/data/fastqLinks.csv"
    output: temp(workDir + "/Results/temp_fastq/runs/{name}.fq.gz")
    params:
        partition = getPartition,
        link = getLink
    resources:
        mem_mb = 3000,
        cpus = 1,
        time = 120
    shell:
     """
     wget -O {output} {params.link} 
     """
     
rule mergeFASTQ:
    input:
        r1 = lambda wildcards: expand(workDir + "/Results/temp_fastq/runs/{run}_1.fq.gz", run = SAMPLES_RUN_DICT[wildcards.sample]),
        r2 = lambda wildcards: expand(workDir + "/Results/temp_fastq/runs/{run}_2.fq.gz", run = SAMPLES_RUN_DICT[wildcards.sample])
    output:
        r1 = workDir + "/Results/temp_fastq/{sample}_R1.fq",
        r2 = workDir + "/Results/temp_fastq/{sample}_R2.fq"
    params: partition = getPartition,
    resources:
        mem_mb = 3000,
        cpus = 1,
        time = 120
    shell:
     """
     zcat {input.r1} > {output.r1}
     zcat {input.r2} > {output.r2}
     """
