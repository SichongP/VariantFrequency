

rule HaplotypeCaller:
    input:
        bam = workDir + "/Results/bams/{sample}.markDup.sorted.bam",
        fa = workDir + "/data/regionsRef.fa"
    output: workDir + "/Results/GVCFs/{sample}.g.vcf.gz"
    params: tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    shell:
     """
     {params.tool} --java-options "-Xmx4g" HaplotypeCaller -R {input.fa} -I {input.bam} -O {output} -ERC GVCF
     """

rule combineGVCF:
    input:
        gvcf = expand(workDir + "/Results/GVCFs/{sample}.g.vcf.gz", sample = SAMPLES) ,
        fa = workDir + "/data/regionsRef.fa"
    output: workDir + "/Results/GVCFs/combined.g.vcf.gz"
    params: tool = workDir + "/tools/gatk-4.2.0.0/gatk", inputGVCF = lambda wildcards, input: ["-V " + file for file in input['gvcf']]
    shell:
     """
     {params.tool} CombineGVCFs -R {input.fa} {params.inputGVCF} -O {output}
     """

rule genotypeGVCF:
    input:
        vcf = workDir + "/Results/GVCFs/combined.g.vcf.gz",
        fa = workDir + "/data/regionsRef.fa",
        targetVCF = workDir + "/data/remappedVariants.vcf"
    output: workDir + "/Results/VCF/output.vcf.gz"
    params: tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    shell:
     """
     {params.tool} --java-options "-Xmx4g" GenotypeGVCFs -R {input.fa} -V {input.vcf} -L {input.targetVCF} --include-non-variant-sites -O {output}
     """
