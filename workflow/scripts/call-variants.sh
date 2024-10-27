
# Output for GATK preprocessing
OUTPUT="output/gatk"

# Binary alignmnt file.
BAM=${snakemake_output[sorted_bam]}/mapped_reads.sorted.bam

# Processed alignment files.
RG=${snakemake_output[rg]}/mapped_reads.rg.bam
DEDUP=${snakemake_output[dedup]}/mapped_reads.dedup.bam

# Reference FASTA
REF=${snakemake_input["ref"]}

# Variant calling file.
VCF=${OUTPUT}/snps.vcf

# Annotation database.
DB=Geobacter_sulfurreducens_kn400

# SNP calling.
gatk HaplotypeCaller -R ${REF} -I ${DEDUP} -O ${VCF} \
  --sequence-dictionary ${OUTPUT}/ref2.dict

# Annotate variants.
snpEff ${DB} ${VCF} >${OUTPUT}/snps.ann.vcf