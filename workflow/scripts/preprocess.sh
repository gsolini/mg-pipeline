#!/usr/bin/env bash
#
# mg-pipeline module 01 — read preprocessing
#
# Deduplication (optional) -> adapter removal -> quality trimming.
# Standalone: takes everything as arguments, reads no config file.
#
# Usage:
#   preprocess.sh --r1 IN_R1 --r2 IN_R2 --out1 OUT_R1 --out2 OUT_R2 [options]
#
# Required:
#   --r1 FILE          input R1 fastq(.gz)
#   --r2 FILE          input R2 fastq(.gz)
#   --out1 FILE        output R1 fastq.gz
#   --out2 FILE        output R2 fastq.gz
#
# Options (defaults in brackets):
#   --threads N        [8]
#   --adapters PATH    adapter fasta, or "bbmap" for BBTools' bundled set [bbmap]
#   --min-len N        minimum read length after trimming [50]
#   --trimq N          quality trim threshold [20]
#   --k N              kmer length for adapter matching [23]
#   --mink N           minimum kmer at read ends [11]
#   --hdist N          allowed mismatches in adapter kmers [1]
#   --entropy F        low-complexity filter, 0 disables [0.0]
#   --ftm N            force-trim modulo; 5 trims 151bp reads to 150 [5]
#   --dedupe BOOL      true|false, run clumpify before trimming [true]
#   --platform NAME    sequencer, sets optical-duplicate distance [other]
#                        nextseq      dupedist=40, spany adjacent
#                        hiseq        dupedist=40   (1T, 2500)
#                        hiseq_patterned  dupedist=2500  (3000, 4000, X)
#                        novaseq      dupedist=12000
#                        other        PCR duplicates only, no optical pass
#   --tmpdir DIR       scratch for intermediates [system temp]
#   --log FILE         append tool stderr here as well as to this script's stderr
#
# Notes:
#   Deduplication runs BEFORE trimming: clumpify groups reads by sequence
#   similarity, and trimming makes identical fragments look different.
#   For Hi-C reads, pass --dedupe false. Hi-C duplicates must be identified
#   after alignment from the mapped positions of both ends (pairtools dedup);
#   read-level dedup cannot do that, and destroys the duplicate-rate QC metric.

set -euo pipefail

# --- defaults ----------------------------------------------------------------
threads=8
adapters="bbmap"
min_len=50
trimq=20
k=23
mink=11
hdist=1
entropy=0.0
ftm=5
dedupe=true
platform="other"
tmpdir="${TMPDIR:-/tmp}"
logfile=""

# --- parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r1)       r1="$2";       shift 2 ;;
    --r2)       r2="$2";       shift 2 ;;
    --out1)     out1="$2";     shift 2 ;;
    --out2)     out2="$2";     shift 2 ;;
    --threads)  threads="$2";  shift 2 ;;
    --adapters) adapters="$2"; shift 2 ;;
    --min-len)  min_len="$2";  shift 2 ;;
    --trimq)    trimq="$2";    shift 2 ;;
    --k)        k="$2";        shift 2 ;;
    --mink)     mink="$2";     shift 2 ;;
    --hdist)    hdist="$2";    shift 2 ;;
    --entropy)  entropy="$2";  shift 2 ;;
    --ftm)      ftm="$2";      shift 2 ;;
    --dedupe)   dedupe="$2";   shift 2 ;;
    --platform) platform="$2"; shift 2 ;;
    --tmpdir)   tmpdir="$2";   shift 2 ;;
    --log)      logfile="$2";  shift 2 ;;
    -h|--help)  sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; exit 1 ;;
  esac
done

for required in r1 r2 out1 out2; do
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
work=$(mktemp -d "${tmpdir%/}/preprocess.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "[preprocess] input:    $r1 / $r2"
echo "[preprocess] workdir:  $work"
echo "[preprocess] dedupe:   $dedupe (platform: $platform)"

# --- optical duplicate settings by platform ----------------------------------
# dupedist values are from the clumpify.sh documentation.
case "$platform" in
  nextseq)          optical_flags="optical dupedist=40 spany adjacent" ;;
  hiseq)            optical_flags="optical dupedist=40" ;;
  hiseq_patterned)  optical_flags="optical dupedist=2500" ;;
  novaseq)          optical_flags="optical dupedist=12000" ;;
  other)            optical_flags="" ;;
  *) echo "ERROR: unknown --platform: $platform" >&2; exit 1 ;;
esac

# --- adapter reference -------------------------------------------------------
# "adapters" is a BBTools alias for its bundled adapters.fa.
if [[ "$adapters" == "bbmap" ]]; then
  adapter_ref="adapters"
else
  [[ -s "$adapters" ]] || { echo "ERROR: adapter file not found: $adapters" >&2; exit 1; }
  adapter_ref="$adapters"
fi

# --- 1. deduplicate ----------------------------------------------------------
# Runs first: clumpify groups by sequence, and trimming would mask duplicates.
if [[ "$dedupe" == "true" ]]; then
  echo "[preprocess] 1/3 clumpify (dedupe)"
  clumpify.sh \
    in="$r1" in2="$r2" \
    out="$work/dedup_R1.fq.gz" out2="$work/dedup_R2.fq.gz" \
    dedupe $optical_flags \
    threads="$threads"
  step1_r1="$work/dedup_R1.fq.gz"
  step1_r2="$work/dedup_R2.fq.gz"
else
  echo "[preprocess] 1/3 clumpify skipped (--dedupe false)"
  step1_r1="$r1"
  step1_r2="$r2"
fi

# --- 2. adapter removal ------------------------------------------------------
# ktrim=r trims to the right of an adapter match.
# tbo trims by overlap of the read pair; tpe trims both reads to equal length.
# ftm=5 drops the extra base on 151bp reads so lengths are a multiple of 5.
echo "[preprocess] 2/3 bbduk (adapters)"
bbduk.sh \
  in="$step1_r1" in2="$step1_r2" \
  out="$work/adapt_R1.fq.gz" out2="$work/adapt_R2.fq.gz" \
  ref="$adapter_ref" \
  ktrim=r k="$k" mink="$mink" hdist="$hdist" \
  ftm="$ftm" \
  minlen="$min_len" tpe tbo \
  threads="$threads"

# --- 3. quality trimming -----------------------------------------------------
# qtrim=r trims low-quality bases from the right end only;
echo "[preprocess] 3/3 bbduk (quality)"
entropy_flag=""
if [[ "$(echo "$entropy > 0" | bc -l)" -eq 1 ]]; then
  entropy_flag="entropy=$entropy"
fi

bbduk.sh \
  in="$work/adapt_R1.fq.gz" in2="$work/adapt_R2.fq.gz" \
  out="$out1" out2="$out2" \
  qtrim=r trimq="$trimq" minlen="$min_len" \
  $entropy_flag \
  threads="$threads"

# --- sanity check ------------------------------------------------------------
# Exit codes alone do not catch silent truncation.
for f in "$out1" "$out2"; do
  [[ -s "$f" ]] || { echo "ERROR: output is empty: $f" >&2; exit 1; }
done

in_reads=$(( $(gzip -cd "$r1" | wc -l) / 4 ))
out_reads=$(( $(gzip -cd "$out1" | wc -l) / 4 ))
pct=$(awk -v a="$out_reads" -v b="$in_reads" 'BEGIN{printf "%.1f", (b>0 ? 100*a/b : 0)}')
echo "[preprocess] reads: $in_reads -> $out_reads (${pct}% retained)"

if (( out_reads == 0 )); then
  echo "ERROR: all reads were removed" >&2
  exit 1
fi

echo "[preprocess] done"
