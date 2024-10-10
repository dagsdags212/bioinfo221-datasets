#!/usr/bin/env bash

set -eux

# Micromamba environemnt.
ENV=bioinfo221

# Output directory.
OUTPUT=output

# Binary alignmnt file.
BAM=${OUTPUT}/mapped_reads.sorted.bam

# Processed alignment files.
RG=${OUTPUT}/mapped_reads.rg.bam
DEDUP=${OUTPUT}/mapped_reads.dedup.bam

# Variant calling file.
VCF=${OUTPUT}/snps.vcf

# Annotation database.
DB=Geobacter_sulfurreducens_kn400

# Group all reads to a single header ID.
micromamba run -n ${ENV} gatk AddOrReplaceReadGroups INPUT=${BAM} OUTPUT=${RG} \
  SORT_ORDER=coordinate RGID=sam RGLB=lib RGPL=illumina RGSM=Sample1 RGPU=1 CREATE_INDEX=True

# Mark duplicate alignments.
micromamba run -n ${ENV} gatk MarkDuplicates INPUT=${RG} OUTPUT=${DEDUP} \
  METRICS_FILE=${OUTPUT}/mdup.metrics.txt

# Create dictionary for deduplicated file.
if [ ! -f "${OUTPUT}/ref2.dict" ]; then
  micromamba run -n ${ENV} gatk CreateSequenceDictionary REFERENCE=${REF} OUTPUT=${OUTPUT}/ref2.dict
fi

# Reindex file.
micromamba run -n ${ENV} samtools index ${DEDUP}

# SNP calling.
micromamba run -n ${ENV} gatk HaplotypeCaller -R ${REF} -I ${DEDUP} -O ${VCF} \
  --sequence-dictionary ${OUTPUT}/ref2.dict

# Annotate variants.
micromamba run -n ${ENV} snpEff ${DB} ${VCF} >${OUTPUT}/snps.ann.vcf
