#!/usr/bin/env bash

set -eu

source ./verify.sh

# Set dependencies for the script
dependencies=("R")

verify_dependencies ${dependencies}

### DO NOT EDIT THE CODE BEFORE THIS COMMENT ###

# Create directory for GWAS output
mkdir -p output/

announce "Running GWAS analysis"

Rscript ./run_gwas.R
