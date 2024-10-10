#!/usr/bin/env bash

set -uex

# Output directory.
OUTPUT=output
mkdir -p ${OUTPUT}

# Micromamba environment.
ENV=bioinfo221

# Paired-end reads.
R1=reads/reads1.fq
R2=reads/reads2.fq

# Reference file.
REF=reads/ref2.fasta

# Alignment files.
SAM=${OUTPUT}/mapped_reads.sam
BAM=${OUTPUT}/mapped_reads.sorted.bam

# Generate reference index.
micromamba run -n ${ENV} bowtie2-build ${REF} ref
micromamba run -n ${ENV} samtools faidx ${REF}

# Perform read mapping.
micromamba run -n ${ENV} bowtie2 -x ref -1 ${R1} -2 ${R2} -S ${SAM}

# Convert SAM to BAM.
micromamba run -n ${ENV} gatk SortSam INPUT=${SAM} OUTPUT=${BAM} SO=coordinate
