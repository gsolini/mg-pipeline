# mg-pipeline module 01 — read preprocessing
#
# Writes stage 01. Shotgun and Hi-C reads run through the same rule, separated
# by the read_type wildcard; parameters differ via pp(), which falls back to
# the shared preprocess block for anything the hic block does not override.


# read_type takes exactly two values. Without this, the wildcard would also
# match path fragments containing slashes and produce confusing errors.
wildcard_constraints:
    read_type="shotgun|hic",
    sample="[A-Za-z0-9_.-]+",
    read="R1|R2",

rule preprocess:
    """Dedupe (optional) -> adapter removal -> quality trim."""
    input:
        unpack(get_reads),
    output:
        r1=f"{OUTDIR}/01_qc/{{sample}}_{{read_type}}_clean_R1.fq.gz",
        r2=f"{OUTDIR}/01_qc/{{sample}}_{{read_type}}_clean_R2.fq.gz",
    params:
        script=join(SCRIPTS, "preprocess.sh"),
        adapters=lambda wc: pp("adapters", wc.read_type),
        min_len=lambda wc: pp("min_len", wc.read_type),
        trimq=lambda wc: pp("trimq", wc.read_type),
        k=lambda wc: pp("k", wc.read_type),
        mink=lambda wc: pp("mink", wc.read_type),
        hdist=lambda wc: pp("hdist", wc.read_type),
        entropy=lambda wc: pp("entropy", wc.read_type),
        ftm=lambda wc: pp("ftm", wc.read_type),
        # false for Hi-C: duplicates are called post-alignment by pairtools.
        dedupe=lambda wc: str(pp("dedupe", wc.read_type)).lower(),
        platform=lambda wc: sample_platform(wc.sample, wc.read_type),
    log:
        f"{OUTDIR}/logs/preprocess/{{read_type}}_{{sample}}.log",
    benchmark:
        f"{OUTDIR}/benchmarks/preprocess/{{read_type}}_{{sample}}.tsv"
    threads: nthreads("preprocess")
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bash '{params.script}' "
        "--r1 '{input.r1}' --r2 '{input.r2}' "
        "--out1 '{output.r1}' --out2 '{output.r2}' "
        "--threads {threads} "
        "--adapters '{params.adapters}' "
        "--min-len {params.min_len} "
        "--trimq {params.trimq} "
        "--k {params.k} --mink {params.mink} --hdist {params.hdist} "
        "--entropy {params.entropy} "
        "--ftm {params.ftm} "
        "--dedupe {params.dedupe} "
        "--platform {params.platform} "
        "--tmpdir '{resources.tmpdir}' "
        "&> {log}"


rule fastqc:
    """QC report for one preprocessed fastq.

    Separate from preprocessing so reports run in parallel, and so the same
    rule can be pointed at raw reads for a before/after comparison.
    """
    input:
        f"{OUTDIR}/01_qc/{{sample}}_{{read_type}}_clean_{{read}}.fq.gz",
    output:
        html=f"{OUTDIR}/01_qc/fastqc/{{sample}}_{{read_type}}_clean_{{read}}_fastqc.html",
        zip=f"{OUTDIR}/01_qc/fastqc/{{sample}}_{{read_type}}_clean_{{read}}_fastqc.zip",

    params:
        outdir=lambda wc, output: os.path.dirname(output.html),
    log:
        f"{OUTDIR}/logs/fastqc/{{read_type}}_{{sample}}_{{read}}.log",
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc -o {params.outdir} -t {threads} {input} &> {log}"


rule preprocess_all:
    """Convenience target: run module 01 for every sample and stop.

    snakemake --use-conda --cores 8 preprocess_all
    """
    input:
        expand(
            f"{OUTDIR}/01_qc/fastqc/{{sample}}_{{read_type}}_clean_{{read}}_fastqc.html",
            read_type=["shotgun", "hic"],
            sample=SHOTGUN_SAMPLES,
            read=["R1", "R2"],
        ),
        expand(
            f"{OUTDIR}/01_qc/fastqc/{{sample}}_{{read_type}}_clean_{{read}}_fastqc.zip",
            read_type=["shotgun", "hic"],
            sample=HIC_SAMPLES,
            read=["R1", "R2"],
        ),
