# Module 01 — Read preprocessing

Removes duplicates, adapters, and low-quality bases from raw reads, then runs
FastQC on the result. Runs on both shotgun and Hi-C libraries, with different
duplicate handling for each.

Steps, in order:

1. **Deduplicate** (`clumpify.sh`) — optional, skipped for Hi-C
2. **Adapter removal** (`bbduk.sh`) — also force-trims the junk 151st base
3. **Quality trimming** (`bbduk.sh`)
4. **FastQC** on the cleaned reads

Writes stage `01_qc`. Always runs; there is no toggle to disable it.

---

## Inputs

Read from the sample sheet (`config/samples.tsv`, or whatever `samples:` points at):

| Column | Used for |
|---|---|
| `sample` | output directory and file prefix |
| `r1`, `r2` | shotgun reads (`read_type=shotgun`) |
| `hic_r1`, `hic_r2` | Hi-C reads (`read_type=hic`), blank if absent |

Paths may be absolute, relative to `fastq_dir`, or start with `~`. An absolute
path in the sheet overrides `fastq_dir`.

Samples with blank Hi-C columns are skipped by the Hi-C branch rather than
failing.

---

## Outputs

```
<outdir>/01_qc/{read_type}/{sample}/
├── {sample}_clean_R1.fq.gz
├── {sample}_clean_R2.fq.gz
└── fastqc/
    ├── {sample}_clean_R1_fastqc.html
    ├── {sample}_clean_R1_fastqc.zip
    ├── {sample}_clean_R2_fastqc.html
    └── {sample}_clean_R2_fastqc.zip

<outdir>/logs/preprocess/{read_type}_{sample}.log
<outdir>/logs/fastqc/{read_type}_{sample}_{read}.log
<outdir>/benchmarks/preprocess/{read_type}_{sample}.tsv
```

`{read_type}` is `shotgun` or `hic`. The cleaned fastqs are what module 02
(decontamination) and module 03 (assembly) consume.

The preprocessing log ends with a line reporting input reads, output reads, and
percent retained. Healthy shotgun libraries retain roughly 85–95%. The bundled 
test dataset (biosolids, 10M pairs) retains 95.9%. Much lower on your own data 
usually means `trimq` or `min_len` is too aggressive.

---

## Configuration

All keys live under `preprocess:` in the config.

| Key | Default | Meaning |
|---|---|---|
| `adapters` | `bbmap` | Adapter FASTA, or `bbmap` for BBTools' bundled set |
| `min_len` | `50` | Minimum read length after trimming (bp) |
| `trimq` | `20` | Quality threshold for trimming (phred) |
| `k` | `23` | Kmer length for adapter matching |
| `mink` | `11` | Minimum kmer length at read ends |
| `hdist` | `1` | Allowed mismatches in adapter kmers |
| `entropy` | `0.0` | Low-complexity filter; 0 disables |
| `ftm` | `5` | Force-trim modulo; 5 trims 151bp reads to 150. 0 disables |
| `dedupe` | `true` | Run clumpify before trimming |
| `platform` | `other` | Sequencer, selects optical-duplicate settings |

### Platform values

`platform` sets clumpify's `dupedist`, taken from the clumpify documentation:

| Value | Flags | Instruments |
|---|---|---|
| `nextseq` | `optical dupedist=40 spany adjacent` | NextSeq (tile-edge smearing) |
| `hiseq` | `optical dupedist=40` | HiSeq 1T, HiSeq 2500 |
| `hiseq_patterned` | `optical dupedist=2500` | HiSeq 3000, 4000, X |
| `novaseq` | `optical dupedist=12000` | NovaSeq |
| `other` | *(none)* | PCR duplicates only, no optical pass |

`other` is the default because the optical pass needs flowcell coordinates in
the read headers; applying it where they are absent does nothing useful.

### Hi-C overrides

Keys under `preprocess: hic:` apply only to Hi-C libraries; anything not listed
there is inherited from the shared block. Currently:

```yaml
preprocess:
  hic:
    dedupe: false
```

Deleting a key from the `hic:` block reverts it to the shared value.

---

## Running it

### As part of the pipeline

```bash
# every sample, both read types
snakemake --use-conda --cores 8 preprocess_all

# one specific output
snakemake --use-conda --cores 8 \
  results/01_qc/shotgun/BS_01/BS_01_clean_R1.fq.gz

# on the test dataset
snakemake --configfile test/config.yaml --use-conda --cores 4 preprocess_all
```

### Standalone

`workflow/scripts/preprocess.sh` takes everything as arguments and reads no
config file, so it runs anywhere BBTools is available:

```bash
workflow/scripts/preprocess.sh \
  --r1 raw/sample_R1.fq.gz --r2 raw/sample_R2.fq.gz \
  --out1 clean/sample_R1.fq.gz --out2 clean/sample_R2.fq.gz \
  --threads 8 --platform nextseq
```

`--help` lists every option. For Hi-C libraries, add `--dedupe false`.

---

## Notes and gotchas

**Deduplication runs before trimming.** Clumpify groups reads by sequence
similarity, so trimming first makes identical fragments look different and
duplicates get missed.

**Hi-C reads are never deduplicated here.** Hi-C duplicates must be identified
after alignment, from the mapped positions of both ends (`pairtools dedup`).
Read-level dedup cannot do that correctly, and it destroys the duplicate-rate
QC metric — which for Hi-C is a readout of prep quality, not just noise.
Hi-C libraries routinely run 80–95% duplicates by design.

**`ftm` is applied in the adapter step, not the quality step.** BBTools applies
force-trim operations before kmer matching, so `ftm=5` removes the junk 151st
base from raw reads and leaves adapter-trimmed reads at their natural length.
In the quality step it would instead round every already-trimmed read down to
a multiple of 5, costing up to 4 bases per read.

**Intermediates go to a temp directory that is cleaned up on exit**, including
on failure, via `mktemp -d` and a `trap`. The location comes from Snakemake's
`{resources.tmpdir}`, set in `profiles/slurm/config.yaml` for cluster runs.

**FastQC is a separate rule** so reports run in parallel with other work, and
so the same rule can be pointed at raw reads for a before/after comparison.

---

## Tools

| Tool | Env | Purpose |
|---|---|---|
| BBTools (`clumpify.sh`, `bbduk.sh`) | `workflow/envs/bbmap.yaml` | dedup, trimming |
| FastQC | `workflow/envs/fastqc.yaml` | QC reports |

On Apple Silicon, bioconda may have no `osx-arm64` build. Prefix the conda or
snakemake command with `CONDA_SUBDIR=osx-64` to build under Rosetta.
