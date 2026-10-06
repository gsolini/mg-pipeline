# mg-pipeline — shared setup and helpers
#
# Included by the Snakefile AFTER `configfile:`, and BEFORE any rule file.
# Contains no rules: only the sample sheet and functions the rules call.

import os
from os.path import join, expanduser
import pandas as pd


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
    
def sample_platform(sample):
    """Platform from the sample sheet, falling back to the config default."""
    if "platform" in samples.columns:
        val = samples.loc[sample, "platform"]
        if pd.notna(val) and str(val).strip():
            return str(val).strip()
    return config["preprocess"]["platform"]


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
