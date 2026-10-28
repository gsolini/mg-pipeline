#!/usr/bin/env bash
#
# mg-pipeline module 02 — host decontamination
#
# Standalone: takes everything as arguments, reads no config file.
#
# Usage:
#   decontaminate.sh --r1 IN_R1 --r2 IN_R2 --out1 OUT_R1 --out2 OUT_R2 [options]
#
# Required:
#   --r1 FILE          input R1 fastq(.gz)
#   --r2 FILE          input R2 fastq(.gz)
#   --out1 FILE        output R1 fastq.gz
#   --out2 FILE        output R2 fastq.gz
#	--index FILE       host genome index (.bt2)
#
# Options (defaults in brackets):
#   --threads N        [8]
#   --stats FILE       output file to write hostile's decontam stats to
#   --aligner          {bowtie2,minimap2,auto}
#   --cache DIR        alternate cache directory
#   --tmpdir DIR       scratch for intermediates [system temp]
#   --log FILE         append tool stderr here as well as to this script's stderr
#

set -euo pipefail

# --- defaults ----------------------------------------------------------------
threads=8
aligner="auto"
tmpdir="${TMPDIR:-/tmp}"
logfile=""
cache=""
stats=""

# --- parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r1)   	r1="$2";       	shift 2 ;;
    --r2)   	r2="$2";       	shift 2 ;;
    --out1)     out1="$2";     	shift 2 ;;
    --out2)     out2="$2";     	shift 2 ;;
	--index)	index="$2";    	shift 2 ;;
    --threads)  threads="$2";  	shift 2 ;;
    --aligner) 	aligner="$2";  	shift 2 ;;
    --tmpdir)   tmpdir="$2";   	shift 2 ;;
    --log)      logfile="$2";  	shift 2 ;;
	--cache)	cache="$2";  	shift 2 ;;
	--stats)    stats="$2";    shift 2 ;;
    -h|--help)  sed -n '2,50p' "$0"; exit 0 ;;					## update to reflect line numbers of new script
    *) echo "ERROR: unknown argument: $1" >&2; exit 1 ;;
  esac
done

for required in r1 r2 out1 out2 index; do
  if [[ -z "${!required:-}" ]]; then
    echo "ERROR: --${required} is required" >&2
    exit 1
  fi
done

for f in "$r1" "$r2"; do
  [[ -s "$f" ]] || { echo "ERROR: input not found or empty: $f" >&2; exit 1; }
done

if [[ -n "$logfile" ]]; then
  mkdir -p "$(dirname "$logfile")"
  exec 2> >(tee -a "$logfile" >&2)
fi

mkdir -p "$(dirname "$out1")" "$(dirname "$out2")" "$tmpdir"

# --- working directory, cleaned up on any exit -------------------------------
work=$(mktemp -d "${tmpdir%/}/decontam.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "[decontaminate] input:    $r1 / $r2"
echo "[decontaminate] workdir:  $work"
echo "[decontaminate] index:    $index"

# --- 1. Decontaminate with hostile -------------------------------------------
# Runs first: depending upon whether you are using hostile's built-in index or
# providing your own custom genome, hostile will either download the pre-built
# index or create one from your custom input
if [[ -n "$cache" ]]; then
  export HOSTILE_CACHE_DIR="$cache"
fi

hostile clean \
  --fastq1 "$r1" --fastq2 "$r2" \
  --index "$index" \
  --output "$work/hostile" \
  --threads "$threads" \
  --aligner "$aligner" \
  > "$work/hostile.json"

# Hostile reports its own output paths; do not infer them from the input names.
cleaned1=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))[0]['fastq1_out_path'])" "$work/hostile.json")
cleaned2=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))[0]['fastq2_out_path'])" "$work/hostile.json")

mv "$cleaned1" "$out1"
mv "$cleaned2" "$out2"

# Write out the stats file if specified
if [[ -n "$stats" ]]; then
  mkdir -p "$(dirname "$stats")"
  cp "$work/hostile.json" "$stats"
fi

# --- sanity check ------------------------------------------------------------
# Exit codes alone do not catch silent truncation.
for f in "$out1" "$out2"; do
  [[ -s "$f" ]] || { echo "ERROR: output is empty: $f" >&2; exit 1; }
done

removed=$(python -c "import json,sys; print(f\"{json.load(open(sys.argv[1]))[0]['reads_removed_proportion']:.2%}\")" "$work/hostile.json")
echo "[decontam] host reads removed: $removed"
echo "[decontam] done"


