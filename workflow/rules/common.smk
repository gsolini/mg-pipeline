# mg-pipeline — shared setup and helpers
#
# Included by the Snakefile AFTER `configfile:`, and BEFORE any rule file.
# Contains no rules: only the sample sheet and functions the rules call.

import os
from os.path import join, expanduser
import pandas as pd
import re


# === Paths ===================================================================

# Directory holding workflow/scripts, resolved from the Snakefile's location so
# the pipeline works regardless of the directory Snakemake was launched from.
SCRIPTS = join(workflow.basedir, "scripts")

OUTDIR = config["outdir"]


def scratch(*parts):
    """Path under scratch_dir, fall back to <outdir>/intermediates.

    scratch_dir holds large regenerable intermediates. An empty value in the
    config means "not set", not "the filesystem root".
    """
    base = config.get("scratch_dir") or join(OUTDIR, "intermediates")
    return join(base, *parts)


def resolve(path):
    """Resolve one sample-sheet path against fastq_dir.

    Accepts absolute paths (which win over fastq_dir), paths relative to
    fastq_dir, and paths beginning with ~.
    """
    return join(expanduser(config.get("fastq_dir", "")), expanduser(str(path)))


# === Sample sheet ============================================================

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
samples = samples.apply(lambda c: c.str.strip() if c.dtype == "object" else c)

def _validate_samples(df):
    """Check the sample sheet before any rule runs."""
    required = ["sample", "r1", "r2"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"sample sheet is missing required column(s): {missing}")

    dupes = df["sample"][df["sample"].duplicated()].tolist()
    if dupes:
        raise ValueError(f"duplicate sample IDs in sample sheet: {sorted(set(dupes))}")

    for s in df["sample"]:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(s)):
            raise ValueError(
                f"sample ID '{s}' contains characters that break path wildcards; "
                "use letters, digits, underscore, dot or hyphen only"
            )

    # Hi-C must be paired or absent, never half-specified.
    if "hic_r1" in df.columns or "hic_r2" in df.columns:
        for _, row in df.iterrows():
            s = row["sample"]
            if str(row["r1"]).strip() == str(row["r2"]).strip():
                raise ValueError(f"sample {s}: r1 and r2 are the same file")
            one = str(row.get("hic_r1", "") or "").strip()
            two = str(row.get("hic_r2", "") or "").strip()
            if bool(one) != bool(two):
                raise ValueError(f"sample {s}: hic_r1 and hic_r2 must both be set or both blank")
            if one and one == two:
                raise ValueError(f"sample {s}: hic_r1 and hic_r2 are the same file")

_validate_samples(samples)

# Optional subset filter, for running one or a few samples:
#   snakemake --config only_samples=BS_01 preprocess_all
_only = config.get("only_samples")
if _only:
    keep = {s.strip() for s in str(_only).split(",")}
    missing = keep - set(samples.index)
    if missing:
        raise ValueError(f"only_samples: not in sample sheet: {sorted(missing)}")
    samples = samples.loc[sorted(keep)]

# Samples whose sheet row carries Hi-C reads.
HIC_SAMPLES = [
    s for s in samples.index
    if pd.notna(samples.loc[s, "hic_r1"]) and samples.loc[s, "hic_r1"].strip()
]

SHOTGUN_SAMPLES = list(samples.index)

# group -> [samples]. Samples with no group are binned alone and are absent here.
GROUPS = {
    g: sorted(rows.index)
    for g, rows in samples.dropna(subset=["group"]).groupby("group")
}


def sample_group(sample):
    """The sample's group, or None when it is ungrouped."""
    g = samples.loc[sample, "group"]
    return g if pd.notna(g) and str(g).strip() else None

def sample_platform(sample, read_type="shotgun"):
    """Platform for one library.

    Resolution order, first hit wins:
      1. hic_platform column (Hi-C only) — for the rare split-instrument case
      2. platform column
      3. preprocess.hic.platform in the config (Hi-C only)
      4. preprocess.platform in the config
    """
    def col(name):
        if name in samples.columns:
            val = samples.loc[sample, name]
            if pd.notna(val) and str(val).strip():
                return str(val).strip()
        return None

    if read_type == "hic":
        return col("hic_platform") or col("platform") or pp("platform", "hic")
    return col("platform") or pp("platform", "shotgun")

# === Input functions =========================================================

def get_reads(wildcards):
    """Raw R1/R2 for one sample, picking the shotgun or Hi-C columns.

    Used with unpack(), so the rule sees {input.r1} and {input.r2}.
    """
    row = samples.loc[wildcards.sample]
    if wildcards.read_type == "hic":
        r1, r2 = row["hic_r1"], row["hic_r2"]
        if pd.isna(r1) or not str(r1).strip():
            raise ValueError(f"sample {wildcards.sample} has no Hi-C reads in the sample sheet")
    else:
        r1, r2 = row["r1"], row["r2"]
    return {"r1": resolve(r1), "r2": resolve(r2)}


# === Config accessors ========================================================

def pp(key, read_type):
    """A preprocess param, with Hi-C overrides falling back to shared defaults.

    config["preprocess"]["hic"] lists only the keys that differ for Hi-C;
    anything absent there is inherited from the shared block.
    """
    block = config["preprocess"]
    if read_type == "hic":
        return block.get("hic", {}).get(key, block[key])
    return block[key]


def nthreads(rule_name):
    """Thread count for a rule, from the config's threads block."""
    return config["threads"][rule_name]
