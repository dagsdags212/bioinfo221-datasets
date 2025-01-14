#!/usr/bin/env bash

set -eu

source ./verify.sh

usage() {
  echo ""
  echo "usage: ./align.sh <INPUT_FA> <OUTPUT_PATH>"
  echo ""
}

# Set dependencies for the script
dependencies=("mafft" "raxmlHPC" "jalview")

verify_dependencies ${dependencies}

### DO NOT EDIT THE CODE BEFORE THIS COMMENT ###

# Path to input alignment file
FA=data/ncov.fasta

# Path to output file
OUT=ncov.aln.fasta

# Set MAFFT options
MAXITER=1000
GLOBAL=true

# Set RAxML parameters
PREFIX=$(basename -s .fasta ${OUT})
SEED=12345
MODEL=GTRCAT
THREADS=2

run_mafft() {
  # Create output directory
  mkdir -p output/mafft
  # Choose alignment type
  if [ "${GLOBAL}" == "true" ]; then
    mafft_opts="--globalpair --maxiterate ${MAXITER}"
  else
    mafft_opts="--local --maxiterate ${MAXITER}"
  fi
  # Compose MAFFT command
  # mafft_cmd="mafft ${mafft_opts} ${FA} > output/mafft/${OUT}"
  mafft_cmd="mafft ${FA} > output/mafft/${OUT}"
  # Run MAFFT
  eval ${mafft_cmd}
}

# Generate an HTML file for visualizing MAFFT output
render_msa() {
  jalview --open output/mafft/${OUT} --color clustal --format html --image output/${PREFIX}.html --headless --close
}

run_raxml() {
  # Create output directory
  mkdir -p output/raxml
  # compose RAxML command
  raxml_cmd="raxmlHPC -s output/mafft/${OUT} -n ${PREFIX} -p ${SEED} -m ${MODEL} -T ${THREADS}"
  # Run RAxML
  eval ${raxml_cmd}
  mv RAxML* output/raxml/
}

# Run MAFFT alignment
run_mafft

# Visualize MAFFT results
render_msa

# Ru RAxML tree inference
run_raxml
