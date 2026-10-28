# Module 02 — Host decontamination

Maps preprocessed shotgun reads against a masked host index with
[Hostile](https://github.com/bede/hostile) and discards anything that aligns.
Optional, and off by default.

Steps:

1. **Fetch or build the host index** — once per run, not per sample
2. **Align and filter** (`hostile clean`) — Bowtie2 under the hood for short reads
3. **Record what was removed** — a JSON stats file per sample

Writes stage `02_host_removed`. Shotgun only; see the gotchas below.

---

## Inputs

Preprocessed shotgun reads from stage 01:

```
<outdir>/01_qc/{sample}_shotgun_clean_R1.fq.gz
<outdir>/01_qc/{sample}_shotgun_clean_R2.fq.gz
```

Plus the host index, fetched or built by `decontam_fetch_index` /
`decontam_build_index` into `<resources_dir>/hostile/`.

---

## Outputs

```
<outdir>/02_host_removed/
├── {sample}_shotgun_decon_R1.fq.gz
├── {sample}_shotgun_decon_R2.fq.gz
└── {sample}_shotgun_hostile.json

<outdir>/logs/decontam/shotgun_{sample}.log
<outdir>/benchmarks/decontam/shotgun_{sample}.tsv

<resources_dir>/hostile/{index}/.fetched      # standard index
<resources_dir>/hostile/custom/custom.*.bt2   # custom index
```

The JSON is Hostile's own report. The fields worth reading are `reads_in`,
`reads_out`, `reads_removed` and `reads_removed_proportion` — this is the
record of how much data the stage deleted, so it is kept rather than discarded.

The log also ends with a one-line summary of the proportion removed.

---

## Configuration

Under `decontam:` in the config.

| Key | Default | Meaning |
|---|---|---|
| `enabled` | `false` | Whether the stage runs at all |
| `index` | `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` | Standard Hostile index name |
| `custom_fasta` | `""` | Path to a FASTA; builds a Bowtie2 index from it instead |
| `aligner` | `auto` | `auto` picks Bowtie2 for short reads |

### Standard indexes

Run `hostile index list` for the authoritative set — the list in `config.yaml`
is a convenience copy and can go stale. The ones most likely to matter:

| Index | Masked against |
|---|---|
| `human-t2t-hla` | nothing |
| `human-t2t-hla-argos985` | 985 FDA-ARGOS bacterial genomes |
| `human-t2t-hla.rs-viral-202401_ml-phage-202401` | RefSeq viral + Millard Lab phage |
| `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` | bacteria, viruses and phage |
| `mouse-mm39` | nothing |
| `test-human-mit` | human mitochondrion only; for testing |

Indexes download on first use and are cached under `resources_dir`, so the
download happens once per machine rather than once per run.

### Custom references

```yaml
decontam:
  enabled: true
  custom_fasta: /path/to/host_genome.fa
```

A Bowtie2 index is built from the FASTA and used in place of a standard index.
**The result is unmasked** — run `hostile mask` yourself first if that matters
for your host, and point `custom_fasta` at the masked output.

---

## Running it

```bash
# every shotgun sample
snakemake --use-conda --cores 8 decontam_all

# on the test dataset
snakemake --configfile test/config.yaml --use-conda --cores 4 decontam_all

# one sample
snakemake --config only_samples=ZY_0 --use-conda --cores 8 decontam_all

# check the plan first
snakemake -n decontam_all
```

`decontam_all` runs the stage even when `enabled: false` — naming the target
is treated as an explicit request. The toggle governs whether downstream
modules read from stage 02 or fall back to stage 01.

### Standalone

```bash
workflow/scripts/decontam.sh \
  --r1 clean_R1.fq.gz --r2 clean_R2.fq.gz \
  --out1 decon_R1.fq.gz --out2 decon_R2.fq.gz \
  --index human-t2t-hla-argos985 \
  --threads 8
```

`--help` lists every option. Without `--cache`, Hostile uses its own default
cache location rather than `resources_dir`.

---

## Notes and gotchas

**Hi-C reads are never decontaminated.** Hi-C libraries contain chimeric reads
spanning a ligation junction, and a read whose two halves come from host and
microbe is handled inconsistently by read-level host filters — no current tool
does this correctly. Host DNA in Hi-C is handled after binning instead. This is
a deliberate omission, not an oversight.

**Masking is the reason to use Hostile rather than rolling your own.** An
unmasked human reference discards microbial reads that match conserved human
regions, silently. The default index here is masked against bacteria, viruses
and phage; the phage masking matters specifically because geNomad calls phage
downstream, and losing those reads would be invisible.

**T2T, not hg38.** The default indexes are built on T2T-CHM13v2.0. The regions
absent from GRCh38 are highly repetitive, and reads from them survive host
removal against hg38 and reappear downstream looking like novel microbial
sequence.

**The output paths in the JSON are stale.** Hostile reports where *it* wrote
the files, which is inside a temp directory the script then deletes after
moving the outputs into place. The read counts are accurate; the paths are not.

**The index is built once, not per sample.** It is a separate rule whose output
every decontamination job depends on, which avoids four jobs racing to build
the same index.

**Nothing removed is not the same as nothing happened.** On a sample with no
host DNA — a mock community, say — `reads_removed_proportion` is 0.0 and that
is correct. Check the JSON exists and the index name is the one you expected
before concluding the stage did its job.

---

## Tools

| Tool | Env | Purpose |
|---|---|---|
| Hostile | `workflow/envs/hostile.yaml` | host read removal |
| Bowtie2 | (bundled with Hostile) | alignment, custom index building |

Hostile is a Python package, so the env supplies the Python used to parse its
JSON output.
