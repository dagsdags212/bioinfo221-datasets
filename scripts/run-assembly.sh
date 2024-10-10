#!/usr/env/bin bash

set -eux

# FASTA files.
R1=reads1.fq
R2=reads2.fq
REF=ref2.fasta

# Micromamba environment.
ENV=bioinfo221

# Output file for assembled contigs.
CONTIGS=assembly.out/scaffolds.fasta

# Run assembly.
micromamba run -n ${ENV} spades.py -1 ${R1} -2 ${R2} -o assembly.out

# Compute contigs statistics to assess assembly quality.
micromamba run -n ${ENV} quast -r ${REF} ${CONTIGS}
