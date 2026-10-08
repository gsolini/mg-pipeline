# The sample sheet

The sample sheet is the one file you must write for your own data. It is a
tab-separated table with one row per sample, listing where the reads live and
how samples relate to each other. Everything else the pipeline needs comes from
the config.

A minimal sheet is three columns:

```
sample	r1	r2
ZY_0	ZY_0_R1.fq.gz	ZY_0_R2.fq.gz
ZY_1	ZY_1_R1.fq.gz	ZY_1_R2.fq.gz
```

`config/samples.tsv` in this repository is a template showing every column.
Copy it, delete the columns you do not need, and point `samples:` in your run
config at your copy.

---

## Columns

| Column | Required | Meaning |
|---|---|---|
| `sample` | **yes** | Unique ID. Becomes output directory and file names. |
| `r1` | **yes** | Shotgun forward reads (fastq, usually gzipped). |
| `r2` | **yes** | Shotgun reverse reads. |
| `group` | no | Samples from the same community. Blank = binned alone. |
| `hic_r1` | no | Hi-C forward reads. Blank = no Hi-C for this sample. |
| `hic_r2` | no | Hi-C reverse reads. Must be set if `hic_r1` is. |
| `platform` | no | Sequencer for this sample. Falls back to the config. |
| `hic_platform` | no | Sequencer for the Hi-C reads, when it differs. |

Column order does not matter; columns are matched by their header name. Omit
any optional column entirely and the pipeline uses its fallback. Leave a cell
blank (two tabs in a row) to use the fallback for one sample only.

### `sample`

Must be unique, and must contain only letters, digits, underscores, dots and
hyphens. Spaces and slashes break the output paths and are rejected at startup.

This ID appears throughout the output tree, so pick something you will still
recognise later: `ZY_0`, `WW_A1_2024`, `MSH_fumarole_3`.

### `r1` and `r2`

Paths may be:

- **absolute** — `/central/groups/lab/data/sample_R1.fq.gz`
- **relative to `fastq_dir`** — `sample_R1.fq.gz`, with `fastq_dir` set in the config
- **starting with `~`** — expanded to your home directory

An absolute path always wins, even when `fastq_dir` is set, so the two styles
can be mixed in one sheet. If `fastq_dir` is empty (the default), relative
paths resolve against the directory you launch Snakemake from — usually not
what you want, so use absolute paths or set `fastq_dir`.

### `group`

Samples sharing a group value are treated as the same community: they can be
co-assembled, and every sample in the group is mapped against every assembly in
the group for differential-coverage binning.

Leave it blank for samples that should be binned on their own coverage. Samples
in *different* groups are never cross-aligned — this is what keeps unrelated
datasets (soil vs wastewater) from being compared.

```
sample	group	r1	r2
WW_A1	ww_2024	...	...
WW_B1	ww_2024	...	...     ← binned together, differential coverage
SOIL_1	soil_2024	...	...
SOIL_2	soil_2024	...	...     ← binned together, separately from the above
MSH_01		...	...             ← binned alone
```

Grouping drives the n² cross-mapping, so a group of 10 samples means 100
alignment jobs. Group deliberately.

### `hic_r1` and `hic_r2`

Fill these in only for samples with matched Hi-C libraries. Both must be set or
both blank; half-specified rows are rejected at startup.

Samples without Hi-C simply skip the Hi-C rules — no error, no empty outputs.
Hi-C binners (metaCC, bin3C) run only for samples that have these columns.

### `platform` and `hic_platform`

The sequencer determines how optical duplicates are detected. See the module 01
docs for the accepted values (`nextseq`, `hiseq`, `hiseq_patterned`, `novaseq`,
`other`).

Resolution order for a sample's shotgun reads:

1. `platform` column
2. `preprocess.platform` in the config

For its Hi-C reads:

1. `hic_platform` column
2. `platform` column
3. `preprocess.hic.platform` in the config
4. `preprocess.platform` in the config

**Use the config, not these columns, when the whole dataset shares a platform.**
The ZymoBIOMICS test data is a good example: the WGS run is NovaSeq 6000 and all
the Hi-C runs are NextSeq 500, so `test/config.yaml` sets

```yaml
preprocess:
  platform: novaseq
  hic:
    platform: nextseq
```

and the sheet needs no platform columns at all.

**Use the columns when the platform varies between samples** — most often with
data pulled from SRA, where runs come from different instruments.

---

## Worked examples

### Shotgun only, one group

```
sample	group	r1	r2
WW_A1	ww_2024	WW_A1_R1.fq.gz	WW_A1_R2.fq.gz
WW_B1	ww_2024	WW_B1_R1.fq.gz	WW_B1_R2.fq.gz
WW_C1	ww_2024	WW_C1_R1.fq.gz	WW_C1_R2.fq.gz
```

With `fastq_dir` set in the run config. All three are binned together with
differential coverage.

### Matched Hi-C, mixed grouping

```
sample	group	r1	r2	hic_r1	hic_r2
ZY_0	zymo	ZY_0_R1.fq.gz	ZY_0_R2.fq.gz	ZY_hic_0_R1.fq.gz	ZY_hic_0_R2.fq.gz
ZY_1	zymo	ZY_1_R1.fq.gz	ZY_1_R2.fq.gz	ZY_hic_1_R1.fq.gz	ZY_hic_1_R2.fq.gz
CTRL		CTRL_R1.fq.gz	CTRL_R2.fq.gz
```

`ZY_0` and `ZY_1` are binned together with both coverage and Hi-C binners;
`CTRL` is binned alone with coverage binners only.

### Mixed platforms from SRA

```
sample	group	platform	r1	r2
SRR001	public	novaseq	SRR001_1.fastq.gz	SRR001_2.fastq.gz
SRR002	public	nextseq	SRR002_1.fastq.gz	SRR002_2.fastq.gz
SRR003	public		SRR003_1.fastq.gz	SRR003_2.fastq.gz
```

`SRR003` has no platform listed, so it falls back to `preprocess.platform`.

---

## Common mistakes

**Spaces instead of tabs.** The file must be tab-separated. Many editors convert
tabs to spaces by default. Check with:

```bash
cat -vet config/samples.tsv | head -3
```

Tabs show as `^I`. If you see runs of spaces, fix your editor settings.

**Blank cells that are not actually blank.** A cell containing a space is not
empty, and will be read as a path or group name. The validator catches some of
these, not all.

**Trailing empty rows.** A final newline is fine; a line of nothing but tabs
creates a sample with an empty ID.

**Half-specified Hi-C.** Setting `hic_r1` without `hic_r2` is rejected at
startup rather than failing later during alignment.

**Reusing a sample ID across groups.** IDs must be unique across the entire
sheet, not just within a group.

---

## Checking your sheet

A dry run parses and validates the sheet without doing any work:

```bash
snakemake --configfile my_config.yaml -n preprocess_all
```

Errors in the sheet surface immediately, before any rule executes. To see which
samples the pipeline actually picked up, add `-r` for the full job list.

To run a subset while testing:

```bash
snakemake --configfile my_config.yaml --config only_samples=ZY_0 -n preprocess_all
```
