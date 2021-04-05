

rule HaplotypeCaller:
    input:
        bam = workDir + "/Results/bams/{sample}.markDup.sorted.bam",
        fa = workDir + "/data/regionsRef.fa",
        tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    output: workDir + "/Results/GVCFs/{sample}.g.vcf.gz"
    shell:
     """
     {input.tool} --java-options "-Xmx4g" HaplotypeCaller -R {input.fa} -I {input.bam} -O {output} -ERC GVCF
     """

rule combineGVCF:
    input:
        gvcf = expand(workDir + "/Results/GVCFs/{sample}.g.vcf.gz", sample = SAMPLES) ,
        fa = workDir + "/data/regionsRef.fa",
        tool = workDir + "/tools/gatk-4.2.0.0/gatk"
    output: workDir + "/Results/GVCFs/combined.g.vcf.gz"
    params: inputGVCF = lambda wildcards, input: ["-V " + file for file in input['gvcf']]
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
    output: workDir + "/Results/VCF/output.vcf.gz"
    shell:
     """
     {input.tool} --java-options "-Xmx4g" GenotypeGVCFs -R {input.fa} -V {input.vcf} -L {input.targetVCF} --dbsnp {input.targetVCF} --include-non-variant-sites -O {output}
     """
