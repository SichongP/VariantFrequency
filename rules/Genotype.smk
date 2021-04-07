

rule HaplotypeCaller:
    input:
        bam = workDir + "/Results/bams/{sample}.markDup.sorted.bam",
        bai = workDir + "/Results/bams/{sample}.markDup.sorted.bam.bai",
        fa = workDir + "/data/regionsRef.fa",
        fai = workDir + "/data/regionsRef.fa.fai",
        dict = workDir + "/data/regionsRef.dict,"
        tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    output: workDir + "/Results/GVCFs/{sample}.g.vcf.gz"
    params: partition = getPartition
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 120
    shell:
     """
     {input.tool} --java-options "-Xmx3g" HaplotypeCaller -R {input.fa} -I {input.bam} -O {output} -ERC GVCF
     """

rule combineGVCF:
    input:
        gvcf = expand(workDir + "/Results/GVCFs/{sample}.g.vcf.gz", sample = SAMPLES) ,
        fa = workDir + "/data/regionsRef.fa",
        tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    output: workDir + "/Results/GVCFs/combined.g.vcf.gz"
    params:
        inputGVCF = lambda wildcards, input: ["-V " + file for file in input['gvcf']],
        partition = getPartition
    resources:
        mem_mb = 10000,
        cpus = 1,
        cpus_bmm = 1,
        mem_mb_bmm = 10000,
        time_min = 240
    shell:
     """
     {input.tool} CombineGVCFs -R {input.fa} {params.inputGVCF} -O {output}
     """

rule genotypeGVCF:
    input:
        vcf = workDir + "/Results/GVCFs/combined.g.vcf.gz",
        fa = workDir + "/data/regionsRef.fa",
        targetVCF = workDir + "/data/remappedVariants.vcf",
        tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    output: workDir + "/Results/VCF/output.vcf"
    params: partition = getPartition
    resources:
        mem_mb = 4000,
        cpus = 1,
        time_min = 240
    shell:
     """
     {input.tool} --java-options "-Xmx3g" GenotypeGVCFs -R {input.fa} -V {input.vcf} -L {input.targetVCF} --dbsnp {input.targetVCF} --include-non-variant-sites -O {output}.gz
     gunzip {output}.gz
     """
