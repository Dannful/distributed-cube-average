#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 [-m MIN] [-M MAX] [-o OUTPUT] [-n NP]

Discover SMPI semantic breakpoints via bp_search1 and bp_search2.

Options:
  -m MIN      Minimum message size (default: 1)
  -M MAX      Maximum message size (default: 1000000)
  -o OUTPUT   Output file for breakpoints (default: breakpoints)
  -n NP       Number of MPI processes (default: 2)
  -h          Show this help message
EOF
  exit 1
}

MIN=1
MAX=1000000
OUTPUT="breakpoints"
NP=2

while getopts "m:M:o:n:h" opt; do
  case $opt in
    m) MIN="$OPTARG" ;;
    M) MAX="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    n) NP="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

echo "Running bp_search1..."
BP1=$(mpirun -np "$NP" bp_search1 "$MIN" "$MAX" | tail -1)

echo "Running bp_search2..."
BP2=$(mpirun -np "$NP" bp_search2 "$MIN" "$MAX" | tail -1)

echo "$BP1" > "$OUTPUT"
echo "$BP2" >> "$OUTPUT"

echo "Breakpoints written to $OUTPUT:"
cat "$OUTPUT"
