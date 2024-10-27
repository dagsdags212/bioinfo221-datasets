#!/usr/bin/env bash

set -eux

# Directory for GATK preprocessing output
OUTPUT="output/gatk"

# Reference FASTA
REF=${snakemake_input[ref]}

# Binary alignmnt file.
SORTED_BAM=${snakemake_input[sorted_bam]}

# Output name for sequence dictionary
REF_DICT=${snakemake_input[ref_dict]}

# Processed alignment files.
RG_BAM=${snakemake_output[rg]}
DEDUP_BAM=${snakemake_output[dedup]}
DEDUP_BAM_IDX=${snakemake_output[dedup_idx]}

# Read group parameters
RGID=${snakemake_params["RGID"]}
RGLB=${snakemake_params["RGLB"]}
RGPL=${snakemake_params["RGPL"]}
RGSM=${snakemake_params["RGSM"]}
RGPU=${snakemake_params["RGPU"]}

# Group all reads to a single header ID.
gatk AddOrReplaceReadGroups INPUT=${SORTED_BAM} OUTPUT=${RG_BAM} \
  RGID=${RGID} RGLB=${RGLB} RGPL=${RGPL} RGSM=${RGSM} RGPU=${RGPU} \
  SORT_ORDER=coordinate CREATE_INDEX=True

# Mark duplicate alignments.
gatk MarkDuplicates INPUT=${RG_BAM} OUTPUT=${DEDUP_BAM} \
  METRICS_FILE=${OUTPUT}/mdup.metrics.txt

# Reindex file.
samtools index ${DEDUP_BAM} > ${DEDUP_BAM_IDX}

