# mg-pipeline module 02 — host decontamination
#
# Writes stage 02. Maps preprocessed shotgun reads against a masked host index
# with Hostile and discards anything that aligns.
#
# Shotgun only — Hi-C reads are deliberately skipped, because chimeric reads
# spanning a ligation junction are not handled correctly by read-level host
# filters (see docs/modules/02_decontam.md).
#
# The host index is fetched (standard Hostile index) or built from a custom
# FASTA, once per run, by the two index rules below.

# match path fragments containing slashes and produce confusing errors.
wildcard_constraints:
    index="[A-Za-z0-9_.-]+",

rule decontam_fetch_index:
    """Download a standard Hostile index."""
    output:
        marker=f"{RESOURCES}/hostile/{{index}}/.fetched",
    params:
        cache=f"{RESOURCES}/hostile",
    conda: "../envs/hostile.yaml"
    shell:
        "HOSTILE_CACHE_DIR='{params.cache}' hostile index fetch -n {wildcards.index} "
        "&& touch '{output.marker}'"

rule decontam_build_index:
    """Build a Bowtie2 index from a custom reference FASTA."""
    input:
        fasta=config["decontam"]["custom_fasta"] or [],
    output:
        f"{RESOURCES}/hostile/custom/custom.1.bt2",
    params:
        prefix=f"{RESOURCES}/hostile/custom/custom",
    threads: nthreads("decontam")
    conda: "../envs/hostile.yaml"
    shell:
        "bowtie2-build --threads {threads} '{input.fasta}' '{params.prefix}'"

rule decontam:
    """Remove host reads from shotgun libraries (Hostile).

    Shotgun only. Hi-C reads are deliberately not decontaminated (see
    docs/modules/02_decontam.md).
    """
    input:
        r1=f"{OUTDIR}/01_qc/{{sample}}_shotgun_clean_R1.fq.gz",
        r2=f"{OUTDIR}/01_qc/{{sample}}_shotgun_clean_R2.fq.gz",
        index=(decontam_index_marker()),
    output:
        r1=f"{OUTDIR}/02_host_removed/{{sample}}_shotgun_decon_R1.fq.gz",
        r2=f"{OUTDIR}/02_host_removed/{{sample}}_shotgun_decon_R2.fq.gz",
        stats=f"{OUTDIR}/02_host_removed/{{sample}}_shotgun_hostile.json",
    params:
        script=join(SCRIPTS, "decontam.sh"),
        index=decontam_index_arg(),
        cache=f"{RESOURCES}/hostile",
        aligner=config["decontam"]["aligner"],
    log:
        f"{OUTDIR}/logs/decontam/shotgun_{{sample}}.log",
    benchmark:
        f"{OUTDIR}/benchmarks/decontam/shotgun_{{sample}}.tsv"
    threads: nthreads("decontam")
    conda:
        "../envs/hostile.yaml"
    shell:
        "bash '{params.script}' "
        "--r1 '{input.r1}' --r2 '{input.r2}' "
        "--out1 '{output.r1}' --out2 '{output.r2}' "
        "--index '{params.index}' --cache '{params.cache}' "
        "--aligner '{params.aligner}' "
        "--threads {threads} "
        "--stats '{output.stats}' "
        "&> '{log}'"


rule decontam_all:
    """Convenience target: run module 02 for every shotgun sample and stop.

    snakemake --use-conda --cores 8 decontam_all
    """
    input:
        expand(
            f"{OUTDIR}/02_host_removed/{{sample}}_shotgun_decon_{{read}}.fq.gz",
            sample=SHOTGUN_SAMPLES,
            read=["R1", "R2"],
        ),
        expand(
            f"{OUTDIR}/02_host_removed/{{sample}}_shotgun_hostile.json",
            sample=SHOTGUN_SAMPLES,
        ),



