import os
from glob import glob

configfile: "config.yaml"

INDIR = config["input_reads_dir"]
OUTDIR = config["output_dir"]
LOGDIR = f"{OUTDIR}/logs"
USEARCH_ID = str(int(config["usearch_id"] * 100))

def handle_equiv_ext(path):
    assert os.path.exists(path), f"Input directory does not exist."
    read_abspaths = glob(f"{INDIR}/*.fastq") + glob(f"{INDIR}/*.fq")
    samples, ext = zip(*[(i[0], i[1]) for i in [os.path.splitext(x) for x in read_abspaths]])
    ext = list(set(ext))
    assert len(ext) == 1, f"Please use only .fastq or .fq as the extension of all input read files."
    samples = [os.path.basename(x) for x in samples]
    print(f"{len(samples)} input read files ({ext[0]}) found.\n")
    return samples, ext[0]

SAMPLES, READS_EXT = handle_equiv_ext(INDIR)

rule all:
    input:
        expand(f"{OUTDIR}/{{sample}}/draft/{{sample}}.draft.fa", sample=SAMPLES),
        expand(f"{OUTDIR}/{{sample}}/draft/{{sample}}.usearch_centroids.id{USEARCH_ID}.fa", sample=SAMPLES),
        expand(f"{OUTDIR}/{{sample}}/polish/{{sample}}.final_consensus.fa", sample=SAMPLES),
        f"{OUTDIR}/all_samples.final_consensus.fa"


# runs usearch cluster to generate centroid sequences (representative sequences), then extract the centroid of cluster with highest number of sequences.
rule usearch_centroids_cluster:
    input:
        reads = f"{INDIR}/{{sample}}{READS_EXT}"
    output:
        centroids = f"{OUTDIR}/{{sample}}/draft/{{sample}}.usearch_centroids.id{USEARCH_ID}.fa",
        draft = f"{OUTDIR}/{{sample}}/draft/{{sample}}.draft.fa"
    log:
        f"{LOGDIR}/{{sample}}.usearch_centroids.log"
    shell:
        "export OMP_NUM_THREADS=10; "
        " {config[usearch_path]} -cluster_fast {input.reads} -id {config[usearch_id]} -strand both -centroids {output.centroids} -sizeout 2>&1 | tee {log} "
        ' && awk "/^>/ {{n++}} n>1 {{exit}} {{print}}" {output.centroids} > {output.draft} '


# run medaka polish pipeline
rule medaka_polish_pipeline:
    input:
        reads = f"{INDIR}/{{sample}}{READS_EXT}",
        draft = f"{OUTDIR}/{{sample}}/draft/{{sample}}.draft.fa"
    output:
        polish_outdir = directory(f"{OUTDIR}/{{sample}}/polish"),
        final_consensus = f"{OUTDIR}/{{sample}}/polish/{{sample}}.final_consensus.fa"
    threads: 
        config["threads"]
    log:
        f"{LOGDIR}/{{sample}}.medaka_polish.log"
    conda:
        "envs/medaka.yaml"
    shell:
        "./medaka_polish.sh -r {input.reads} -d {input.draft} -o {output.polish_outdir} -t {threads} -m {config[medaka_model]} 2>&1 | tee {log} "


# merge consensus sequences of all sample into a single FASTA
rule merge_all_consensus:
    input:
        all_consensus = expand(f"{OUTDIR}/{{sample}}/polish/{{sample}}.final_consensus.fa", sample=SAMPLES)
    output:
        merged = f"{OUTDIR}/all_samples.final_consensus.fa"
    params:
        f"{OUTDIR}/*/polish/*.final_consensus.fa"
    shell:
        "cat {params} > {output.merged}"
