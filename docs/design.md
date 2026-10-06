# mg-pipeline — design document

Status: draft, v0.1 (2026-09-28)
Owner: Grace Solini
Purpose: modular, config-driven metagenomics pipeline (preprocess → assembly → binning →
refinement → annotation), version-controlled on GitHub, runnable locally or on Caltech's
Resnick HPC. Successor to the per-project MSH scripts and the pre-model pipeline for metaHiC.

---

## 1. Design goals

1. One config file + one sample sheet drive an entire run, start to finish.
2. Every module is a portable script (bash/python) that can also be run standalone by
   someone who doesn't use the workflow manager.
3. Optional steps (decontamination, Kraken2, Hi-C binning, annotation add-ons) are toggled
   in config, not by editing code.
4. Grouped samples get differential-coverage binning; ungrouped samples get single-sample
   binning.
5. HPC execution is a configuration concern (a SLURM profile), not a separate set of scripts.
6. README/docs are updated as each module is written, not afterwards.

## 2. Framework decision

**Snakemake**, with each module implemented as a standalone script wrapped by a thin rule.

- Local: `snakemake --use-conda --cores 8`
- Resnick: `snakemake --profile profiles/slurm`
- The SLURM profile (partition, memory, time, per-rule resources) replaces a hand-written
  sbatch generator.
- One conda env YAML per tool, so tool versions are pinned and reproducible.

Considered and rejected: hand-chained bash + a custom sbatch generator (reimplements
dependency resolution, resumption, and cluster submission).

Still worth reviewing before heavy investment: **nf-core/mag** (Nextflow) covers a large
portion of steps 1–7 already. The justification for building this pipeline is the Hi-C
binners, geNomad, and the handoff into metaHiC.

## 3. Repository strategy

One repo per tool; analysis projects are separate and *use* the pipeline.

| Repo | Contents |
|---|---|
| `gsolini/metaHiC` (exists) | Hi-C host–MGE model |
| `gsolini/mg-pipeline` (new) | this pipeline |
| per-project analysis repos (e.g. `msh-analysis`) | `config.yaml`, `samples.tsv`, notebooks, figures |

Data and results are never committed. `.gitignore` covers `results/`, `logs/`, `.snakemake/`,
`*.fastq.gz`, `*.bam`, `*.fa`, conda envs, and `.DS_Store`.

## 4. Repository layout

```
mg-pipeline/
├── workflow/
│   ├── Snakefile
│   ├── rules/            # one .smk per module
│   │   ├── preprocess.smk
│   │   ├── decontam.smk
│   │   ├── assembly.smk
│   │   ├── mapping.smk
│   │   ├── binning.smk
│   │   ├── refinement.smk
│   │   └── annotation.smk
│   ├── scripts/          # standalone bash/python; the actual work
│   └── envs/             # one conda YAML per tool
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── profiles/slurm/       # Resnick cluster settings
├── test/                 # tiny subsampled dataset + expected outputs
├── docs/
│   ├── design.md         # this file
│   └── modules/          # one page per module
├── .gitignore
└── README.md
```

## 5. Sample sheet

`config/samples.tsv` — tab-separated. Hi-C columns blank when absent; `group` blank when the
sample is binned alone.

```
sample    group    platform   r1              r2              hic_r1          hic_r2
A1        ww_ucsd   nextseq  A1_R1.fq.gz     A1_R2.fq.gz     A1_hic_R1.fq.gz A1_hic_R2.fq.gz
B1        ww_ucsd   nextseq  B1_R1.fq.gz     B1_R2.fq.gz     B1_hic_R1.fq.gz B1_hic_R2.fq.gz
MSH_01  -   -   MSH_01_R1.fq.gz MSH_01_R2.fq.gz -   -
```
`platform` is optional; samples that leave it blank use preprocess.platform from the config (default is set to `other` unless otherwise specified/changed). (The dashes are for the doc only — in the real TSV those are empty fields between tabs.)

Semantics of `group`: samples in the same group are from the same community and may be
cross-aligned for differential coverage and (optionally) co-assembled. Samples in different
groups must never be cross-aligned.

## 6. Config schema (sketch)

```yaml
samples: config/samples.tsv
outdir: results

preprocess:
  adapters: bbmap
  min_len: 50
  trimq: 20
  k: 23
  mink: 11
  hdist: 1
  entropy: 0.0
  ftm: 5                  # force-trim modulo; applied in the adapter step
  dedupe: true
  platform: other         # fallback when the sample sheet has no platform column

  # Hi-C overrides: only keys that differ from the block above.
  hic:
    dedupe: false         # Hi-C duplicates are called post-alignment by pairtools

decontam:
  enabled: true
  hosts: [human]            # index paths defined in resources

assembly:
  assemblers: [megahit, spades]
  coassembly: false         # per-group co-assembly in addition to per-sample
  selection: n50            # rule module; rules configurable, TBD

kraken2:
  enabled: false

mges:
  genomad: true             # exclude extrachromosomal contigs from binning

binning:
  binners: [metabat2, maxbin2, semibin2, vamb]
  hic_binners: [metacc, bin3c]   # used only when Hi-C columns are present
  differential_coverage: true

refinement:
  binette: true
  drep: true                # across assemblies within a group

annotation:
  gtdbtk: true
  bakta: false
  amrfinder: true
```

## 7. Output tree
```
results/
├── 01_qc/{read_type}/{sample}/                   # read_type = shotgun | hic
├── 02_host_removed/{read_type}/{sample}/
├── 03_assembly/{assembler}/{assembly_id}/        # assembly_id = sample or group_coassembly
├── 04_assembly_qc/{assembly_id}/                 # metaQUAST
│ └── selected/{assembly_id}.fa                   # representative assembly
├── 05_mges/{assembly_id}/                        # geNomad
├── 06a_mapping_shotgun/{assembly_id}/
│ ├── {reads_sample}.bam                          # optional, see mapping_shotgun.keep_bams
│ └── depth.tsv                                   # coverage table consumed by binners
├── 06b_mapping_hic/{assembly_id}/
│ ├── {hic_sample}.bam                            # bwa mem -5SP, chimera-aware
│ └── {hic_sample}.pairs.gz                       # pairtools output; metaCC/bin3C/metaHiC input
├── 07_binning/{assembly_id}/{binner}/
├── 08_refined/{assembly_id}/                     # Binette
├── 09_dereplicated/{group}/                      # dRep
└── 10_annotation/{tool}/{genome}/
logs/{rule}/{wildcards}.log
benchmarks/{rule}/{wildcards}.tsv
```

## 8. Modules

Modules are `.smk` files; stage numbers above are output directories. One module
may write more than one stage.

| Module | File | Writes | Tools | Optional |
|---|---|---|---|
| Read preprocessing | `preprocess.smk` | 01 | bbduk (adapters, quality, length, dedupe) | no |
| Host decontamination | `decontam.smk` | 02 | Hostile or bowtie2 vs host index | yes |
| Assembly + QC + selection | `assembly.smk` | 03, 04 | MEGAHIT, metaSPAdes, metaQUAST, selection | no |
| MGE identification | `mges.smk` | 05 | geNomad | yes |
| Mapping | `mapping.smk` | 06a, 06b | bowtie2/minimap2 (shotgun); bwa + pairtools (Hi-C) | no |
| Binning | `binning.smk` | 07 | MetaBAT2, MaxBin2, SemiBin2, VAMB; metaCC, bin3C | no |
| Refinement | `refinement.smk` | 08, 09 | Binette, dRep | no |
| Annotation | `annotation.smk` | 10 | GTDB-Tk, Bakta, AMRFinderPlus | yes |
| Read Profiling | `profiling.smk` | — | Kraken2 (read-level) | yes |


## 9. Open design decisions

1. **Assembly selection rules.** Which metrics, in what order, with what thresholds. Make
   the selection a script whose rules live in config. Alternative worth considering: bin both
   assemblies and let dRep pick, at the cost of roughly double the binning compute.
2. **Differential-coverage mapping scheme.** Default plan: every sample in a group maps to
   every assembly in that group (n² BAMs per group). Confirm whether VAMB and SemiBin2 should
   instead run in their multi-sample/concatenated-assembly modes, which treat grouping
   differently from MetaBAT2/MaxBin2.
3. **Co-assembly.** Whether per-group co-assembly runs alongside per-sample assemblies by
   default (current wastewater practice is both), and how co-assembly bins merge with
   per-sample bins at the dRep step.
4. **Host decontamination tool.** Replace bmtagger; confirm Hostile vs bowtie2 + index, and
   how custom host indexes (mouse, etc.) are declared in config.
5. **Annotation.** Bakta calls genes internally with Pyrodigal, so it replaces rather than
   supplements Prodigal. Decide whether functional annotation (eggNOG-mapper, DRAM) belongs
   in this pipeline or downstream.
6. **Hi-C handoff.** The mapped Hi-C BAMs (module 05) are the input metaHiC consumes. Define
   that interface explicitly so the two pipelines don't duplicate alignment work.

## 10. Build order

Each step: write → test on `test/` data → document → commit → push.

1. Git setup, repo skeleton, sample-sheet parsing, test dataset (Zymo mock subsampled to
   ~1M reads is a good candidate — known ground truth, small)
2. Preprocessing (01)
3. Decontamination (02)
4. Assembly, QC, selection (03, 04)
5. Mapping (06a, 06b)
6. Binning (07)
7. Binette + dRep (08, 09)
8. Annotation add-ons (10), Kraken2 (11), geNomad (05)

## 11. Config helpers

`workflow/rules/common.smk` holds sample-sheet parsing and shared helpers,
included from the Snakefile after `configfile:`.

Preprocess params fall back from Hi-C overrides to shared defaults:

```python
def pp(key, read_type):
    """Preprocess param, with hic overrides falling back to shared defaults."""
    block = config["preprocess"]
    if read_type == "hic":
        return block.get("hic", {}).get(key, block[key])
    return block[key]
```

Used in rules as `params: trimq=lambda wc: pp("trimq", wc.read_type)`.
The `hic:` block in config.yaml lists only keys that differ from shotgun;
deleting a key from it reverts that param to the shared value.

### Hi-C alignment

Hi-C reads are NOT aligned like shotgun reads. Chimeric reads spanning a
ligation junction are signal, so alignment is single-end and chimera-aware
(`bwa mem -5SP`), then processed by pairtools parse/sort/dedup into a pairs file.

Ported from the existing metaHiC pairtools pipeline — do not rediscover these:
- `pairtools sort --tmpdir` must point at scratch, and `TMPDIR` exported.
  Default `/tmp` on Resnick is RAM-backed and fails on large libraries
  (seen at 404M and 449M reads).
- Exit codes do not catch silent truncation. Validate output size against
  parsed input (sorted output < 1% of input = failure) before the rule's
  output file is written.

Hi-C duplicates are removed here, at the pair level, not during preprocessing.
