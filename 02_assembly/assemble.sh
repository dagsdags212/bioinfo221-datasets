#!/usr/bin/env bash

set -eu

source ./verify.sh

# Set dependencies for the script
dependencies=("spades.py" "quast")

verify_dependencies ${dependencies}

### DO NOT EDIT THE CODE BEFORE THIS COMMENT ###

# Number of cores to use.
THREADS=8

# Path to pair-end FASTQ reads
R1=data/reads1.fq
R2=data/reads2.fq

# Path to reference file
REF=data/ref.fasta

run_spades() {
  mkdir -p output/spades/
  spades_cmd="spades.py -1 ${R1} -2 ${R2} -o output/spades -t ${THREADS}"
  eval ${spades_cmd}
}

run_quast() {
  mkdir -p output/quast/
  quast_cmd="quast output/spades/scaffolds.fasta -o output/quast -r ${REF} -t ${THREADS}"
  eval ${quast_cmd}
}

run_assembly_pipeline() {
  # Run assembly with spades
  run_spades
  # Run assembly evaluation with quast
  run_quast
}

# Run full assembly pipeline
run_assembly_pipeline
