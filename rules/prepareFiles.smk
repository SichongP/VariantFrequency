localrules: getGenomeFASTA, generateRef, generateFASTQLink, indexVCF, indexRef

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
        mapping = workDir + "/data/variantMapping.csv",
        outVCF = workDir + "/data/remappedVariants.vcf"
    conda: workDir + "/envs/python.yaml"
    script: workDir + "/scripts/generateRef.py"

rule indexVCF:
    input: workDir + "/data/remappedVariants.vcf"
    output: workDir + "/data/remappedVariants.vcf.idx"
    params: tool = workDir + "/tools/gatk-4.1.9.0/gatk"
    shell:
     """
     {params.tool} IndexFeatureFile -I data/remappedVariants.vcf
     """

rule indexRef:
    input: workDir + "/data/regionsRef.fa"
    output: workDir + "/data/regionsRef.fa.fai"
    conda: workDir + "/envs/samtools.yaml"
    params: tool = workDir + "/tools/gatk-4.1.9.0/gatk"
    shell:
     """
     samtools faidx {input}
     module load java/1.8
     {params.tool} CreateSequenceDictionary -R {input}
     """
    
rule generateFASTQLink:
    input: accession = workDir + "/data/accessions.tsv"
    output: outLinks = workDir + "/data/fastqLinks.csv"
    conda: workDir + "/envs/python.yaml"
    script: workDir + "/scripts/getFASTQLinks.py"
    
rule getFASTQ:
    output: temp(expand(workDir + "/Results/temp_fastq/runs/{{run}}/{{run}}_{read}.fastq.gz", read = ['1', '2']))
    params:
        partition = getPartition,
        outDir = workDir + "/Results/temp_fastq/runs/"
    resources:
        mem_mb = 3000,
        cpus = 1,
        time_min = 120
    shell:
     """
     bin/enaBrowserTools-1.6/python3/enaDataGet -f fastq -d {params.outDir} {wildcards.run}
     """
     
rule mergeFASTQ:
    input:
        r1 = lambda wildcards: expand(workDir + "/Results/temp_fastq/runs/{run}/{run}_1.fastq.gz", run = SAMPLES_RUN_DICT[wildcards.sample]),
        r2 = lambda wildcards: expand(workDir + "/Results/temp_fastq/runs/{run}/{run}_2.fastq.gz", run = SAMPLES_RUN_DICT[wildcards.sample])
    output:
        r1 = workDir + "/Results/temp_fastq/{sample}_R1.fq",
        r2 = workDir + "/Results/temp_fastq/{sample}_R2.fq"
    params: partition = getPartition,
    resources:
        mem_mb = 3000,
        cpus = 1,
        time_min = 120
    shell:
     """
     zcat {input.r1} > {output.r1}
     zcat {input.r2} > {output.r2}
     """
