#!/usr/bin/env bash

set -eu

source ./verify.sh

# Set dependencies for the script
dependencies=("bowtie2" "samtools" "gatk")

verify_dependencies ${dependencies}

### DO NOT EDIT THE CODE BEFORE THIS COMMENT ###

# Number of cores to use.
THREADS=8

# Path to pair-end FASTQ reads
R1=data/reads1.fq
R2=data/reads2.fq

# Path to reference file
REF=data/ref.fasta

# Prefix for index
PREFIX=ref

# Output filename for alignment file
ALN_OUT=mapped.reads

index_ref() {
  announce "Generating reference index"
  bowtie2-build -f ${REF} data/${PREFIX}
  samtools faidx ${REF}
}

map_reads() {
  announce "Mapping reads to reference"
  # Create directory for bowtie output
  mkdir -p output/bowtie/
  # Generate index for reference file
  bowtie2 -x data/${PREFIX} \
    -1 ${R1} -2 ${R2} \
    -S output/bowtie/${ALN_OUT}.sam
}

sort_sam() {
  announce "Sorting alignment file"
  gatk SortSam \
    INPUT=output/bowtie/${ALN_OUT}.sam OUTPUT=output/gatk/${ALN_OUT}.sorted.bam \
    SO=coordinate
}

add_read_group() {
  announce "Assigning read groups"
  gatk AddOrReplaceReadGroups \
    INPUT=output/gatk/${ALN_OUT}.sorted.bam OUTPUT=output/gatk/${ALN_OUT}.sorted.rg.bam \
    SORT_ORDER=coordinate CREATE_INDEX=True \
    RGID=sam RGLB=lib RGPL=illumina RGSM=Sample1 RGPU=1
}

mark_duplicates() {
  announce "Marking duplicates"
  gatk MarkDuplicates \
    INPUT=output/gatk/${ALN_OUT}.sorted.rg.bam OUTPUT=output/gatk/${ALN_OUT}.sorted.rg.dedup.bam \
    METRICS_FILE=output/gatk/mdup.metrics.txt
}

create_dict() {
  announce "Generating sequence dictionary"
  rm -f output/gatk/ref.dict
  gatk CreateSequenceDictionary \
    REFERENCE=${REF} \
    OUTPUT=output/gatk/ref.dict
  cp output/gatk/ref.dict data/ref.dict
}

index_processed_alignment() {
  announce "Indexing processed BAM file"
  samtools index output/gatk/${ALN_OUT}.sorted.rg.dedup.bam
}

preprocess() {
  announce "Starting data preprocessing"
  mkdir -p output/gatk/
  # Sort mapped reads by genomie coordinate
  sort_sam
  # Assign reads to a read group
  add_read_group
  # Mark mapped duplicates
  mark_duplicates
  # Generate sequence dictionary
  create_dict
  # Index processed BAM file with samtools
  index_processed_alignment
}

call_variants() {
  announce "Calling variants"
  gatk HaplotypeCaller \
    -R ${REF} --sequence-dictionary output/gatk/ref.dict \
    -I output/gatk/${ALN_OUT}.sorted.rg.dedup.bam.bai \
    -O output/gatk/snps.vcf
}

# Generate FM index for reference
#index_ref

# Map reads back to the reference
#map_reads

# Preprocess alignment file
#preprocess

# Identify variants
call_variants
