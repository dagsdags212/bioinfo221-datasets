#!/usr/bin/env bash

set -eu

# Keep track of number of uninstalled programs
errors=0

# Exit with error code 1 if specified program is not installed or in path
check_install() {
  local tool=$1
  if ! command -v ${tool} >/dev/null 2>&1; then
    echo "Error: ${tool} is not installed"
    errors=$((errors + 1))
  fi
}

verify_dependencies() {
  echo "Checking if all dependencies are installed"
  sleep 1
  for dep in ${dependencies[@]}; do
    check_install ${dep}
  done
  if [ "${errors}" -eq 0 ]; then
    echo "All dependencies installed"
  else
    echo "Verify install: ${errors}/${#dependencies[@]} programs was not found in PATH"
    exit 1
  fi
}

announce() {
  echo
  echo "=== $1 ==="
  echo
  sleep 2
}
