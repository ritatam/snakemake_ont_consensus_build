# Consensus sequence construction from ONT long amplicons using USEARCH and medaka - Snakemake workflow
A mini snakemake pipeline for batch constructing consensus sequences from long nanopore amplicons with medaka for a set of samples.

The pipeline includes the following steps:

1. <code>usearch -clutser_fast</code> for clustering read sequences to find centroid with the largest cluster size (i.e. draft sequence)
2. <code>mini_align</code> for aligning the reads to a draft/intermediate sequence
3. <code>medaka consensus</code> for running consensus algorithm across assembly regions
4. <code>medaka stitch</code> for collating results from step 3 to create consensus sequence

Steps 2-4 (wrapped up in <code>medaka_polish.sh</code>) are the medaka polishing steps, and will be iterated for four rounds as advised by medaka documentation. This snakemake script will automatically install medaka package (defined in <code>envs/medaka.yaml</code>) in your working directory without requiring admin priviledges. The final consensus from the 4th polishing round will be renamed as *.final_consensus.fa.

At the end of pipeline, the final consensus sequences of all samples will be merged into a single FASTA file (<code>all_samples.final_consensus.fa</code>) in the output directory.

## Clone this repository

    git clone https://github.com/ritatam/snakemake_ont_consensus_build.git
    cd snakemake_ont_consensus_build

## Dependencies
* conda
* python
* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 
* [medaka](https://github.com/nanoporetech/medaka) 
* [USEARCH](https://www.drive5.com/usearch/download.html) (already part of the release in external_program directory)

## Configuration
Paths and parameters can be edited in <code>config.yaml</code>.
1. <code>input_reads_dir</code>: Input directory containing all the *.fastq amplicon read files (e.g. s_cerevisiae.fastq, a_flavus.fastq, ...) to be processed simultaneously (via snakemake)
2. <code>output_dir</code>: Output directory to store *draft* and *polished* consensus sequences per sample (as subdirectories). All the intermediate files will be stored in the output sample/polish/intermediate_files directory. 
3. <code>usearch_path</code>: Path to the USEARCH binary. Note: This has been already configured to use the USEARCH binary in the repo's external_program folder, so user will not need to download it independently.
4. <code>usearch_id</code>: Minimum sequence identity for USEARCH read clustering, ranging between 0.0 and 1.0.
5. <code>medaka_model</code>: Medaka model for consensus polishing, based on the basecaller used. See medaka [documentation](https://github.com/nanoporetech/medaka#models) for details.
6. <code>threads</code>: Number of threads to use in medaka polishing. 

## Executing the pipeline

Install the required conda environemnt without running the pipeline. Subsequent runs with flag <code>--use-conda</code> will utilise the local environment stored in your working directory (<code>/.snakemake/conda</code>) without requiring internet connection.

    snakemake --use-conda --conda-create-envs-only

Dry-run the pipeline to print the commands and the corresponding i/o files, without running the pipeline.

    snakemake --use-conda -np

Run the pipeline and specify number of cores to use (e.g. 8 cores)

    snakemake --use-conda -c8 
