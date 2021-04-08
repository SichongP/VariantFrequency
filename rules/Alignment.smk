rule indexFASTA:
    input: workDir + "/data/regionsRef.fa"
    output: workDir + "/data/regionsRef.fa.bwt"
    conda: workDir + "/envs/bwa.yaml"
    params: partition = getPartition
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 30
    shell:
     """
     bwa index {input}
     """
rule mapping:
    input: 
        r1 = workDir + "/Results/temp_fastq/{sample}_R1.fq",
        r2 = workDir + "/Results/temp_fastq/{sample}_R2.fq",
        fa = workDir + "/data/regionsRef.fa",
        index = workDir + "/data/regionsRef.fa.bwt"
    output: temp(workDir + "/Results/bams/{sample}.sam")
    resources:
        mem_mb = 8000,
        cpus = 2,
        time_min = 640
    params: partition=getPartition
    conda: workDir + "/envs/bwa.yaml"
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:illumina" -o {output} {input.fa} {input.r1} {input.r2}
     """

rule filterBAM:
    input: workDir + "/Results/bams/{sample}.sam"
    output: temp(workDir + "/Results/bams/{sample}.filtered.bam")
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 240
    params: partition=getPartition
    conda: workDir + "/envs/samtools.yaml"
    shell:
     """
     samtools view -b -F 2820 -o {output} {input}
     """

rule markDuplicate:
    input: workDir + "/Results/bams/{sample}.filtered.bam"
    output: temp(workDir + "/Results/bams/{sample}.markDup.bam")
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 200
    params: partition = getPartition
    conda: workDir + "/envs/sambamba.yaml"
    shell:
     """
     if [ ! -n ${{SLURM_JOBID-}} ]; then
         MYTMPDIR=/scratch/pengsc/$SLURM_JOBID
         cleanup() {{ rm -rf $MYTMPDIR; }}
         trap cleanup EXIT
         mkdir -p $MYTMPDIR
     else
         MYTMPDIR="./temp/"
         mkdir -p $MYTMPDIR
     fi
     sambamba markdup -t {resources.cpus} --tmpdir=$MYTMPDIR {input} {output}
     """

rule sortBAM:
    input: workDir + "/Results/bams/{sample}.markDup.bam"
    output: temp(workDir + "/Results/bams/{sample}.markDup.sorted.bam")
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 240
    params: partition=getPartition
    conda: workDir + "/envs/samtools.yaml"
    shell:
     """
     samtools sort -o {output} {input}
     """

rule indexBAM:
    input: workDir + "/Results/bams/{sample}.markDup.sorted.bam"
    output: temp(workDir + "/Results/bams/{sample}.markDup.sorted.bam.bai")
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 60
    conda: workDir + "/envs/samtools.yaml"
    params: partition = getPartition
    shell:
     """
     samtools index {input}
     """
