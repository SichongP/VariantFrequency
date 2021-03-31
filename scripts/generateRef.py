# This script takes a bed file of variants as input, and generates a ref FASTA file with regions +/- 1kb around the variants, automatically mergeing any regions that overlap
import pandas as pd
import pyranges as pr
from datetime import date

inputBed = snakemake.input['bed']
genomeFa = snakemake.input['ref']
sizeFile = snakemake.input['chromsize']
refFa = snakemake.output['fa']
remappedVar = snakemake.output['remapVar']
vcf = snakemake.output['outVCF']
mapping = snakemake.output['mapping']

try:
    variants = pd.read_csv(inputBed, sep = '\t', header = None)
except:
    print("Error reading variants file: {}. Is it accessible and correct bed format?".format(inputBed))
    exit(1)
chrom_sizes = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in open(sizeFile, 'r').readlines()}
variants.columns = ['Chromosome', 'variantStart', 'variantEnd']

# Expand variant region by 1kb up and down stream
variants['Start'] = variants.apply(lambda row: max(int(row['variantStart']) - 1000, 0), axis = 1)
variants['End'] = variants.apply(lambda row: min(int(row['variantEnd']) + 1000, int(chrom_sizes[row['Chromosome']])), axis = 1)
variants = variants.sort_values(['Chromosome', 'Start', 'End'])

# Cluster overlapping regions
clustered = pr.PyRanges(variants).cluster().df

# Merged overlapping regions
merged = clustered.groupby('Cluster', group_keys = False).apply(lambda df: pd.DataFrame({'Chromosome': df.iloc[0].Chromosome, 'Start': df.Start.min(), 'End': df.End.max()}, index = [df.iloc[0].Cluster]))

# Get new variant positions
clustered['newVariantStart'] = clustered.apply(lambda row: row.variantStart - merged.loc[row.Cluster].Start, axis = 1)
clustered['newVariantEnd'] = clustered.apply(lambda row: row.variantEnd - merged.loc[row.Cluster].Start, axis = 1)
clustered['regionStart'] = clustered.apply(lambda row: merged.loc[row.Cluster].Start, axis = 1)
clustered['regionEnd'] = clustered.apply(lambda row: merged.loc[row.Cluster].End, axis = 1)

# Assign ID
counts = []
for i in range(0,len(clustered)):
    df = clustered[0:i+1]
    counts.append(df.groupby('Cluster')['Chromosome'].transform('size').iloc[-1])
clustered['vid'] = counts
clustered['id'] = 'c' + clustered['Cluster'].astype(str) + 'v' + clustered['vid'].astype(str)

# Get Ref alleles
vcf_df = clustered[['Chromosome', 'variantStart', 'variantEnd']]
vcf_df.columns = ['Chromosome', 'Start', 'End']
ref_alleles = []
for idx, row in vcf_df.iterrows():
    seq = pr.get_fasta(pr.PyRanges(pd.DataFrame(row).transpose()), "/home/pengsc/reference/equcab3/UCSC/equCab3.fa")[0]
    ref_alleles.append(seq.upper())
clustered['ref'] = ref_alleles

# Write to output FASTA:
with open(refFa, 'w') as output:
    for idx, row in merged.iterrows():
        seq = pr.get_fasta(pr.PyRanges(pd.DataFrame(row).transpose()), genomeFa)[0]
        output.write(">{} {}:{}-{}\n".format(str(idx), row.Chromosome, str(row.Start), str(row.End)))
        output.write("{}\n".format(seq))

# Write to output remapped variant list:
clustered[['Chromosome', 'newVariantStart', 'newVariantEnd', 'Cluster', 'id']].to_csv(remappedVar, sep = '\t', index = False, header = False)

#Write to output mapping of old and new variant positions
clustered[['Chromosome', 'variantStart', 'variantEnd', 'newVariantStart', 'newVariantEnd', 'Cluster', 'id']].to_csv(mapping, index = False)

#output to vcf
with open(vcf, 'w') as output:
    output.write('##fileformat=VCFv4.2\n##fileDate={}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'.format(str(date.today().strftime("%Y%m%d"))))
    for idx, row in clustered.iterrows():
        output.write('{chr}\t{pos}\t{id}\t{ref}\t.\t.\t.\t.\n'.format(chr = row.Cluster, pos = str(int(row.newVariantStart) + 1), id = row.id, ref = row.ref))
