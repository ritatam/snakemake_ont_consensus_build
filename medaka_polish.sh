#!/bin/bash 

# This script automates four rounds of medaka polishing workflow.

# The script takes in a nanopore long-read sequencing data file (fastq) and a draft consensus sequence (fasta) to perform four rounds of polishing using the medaka package. The final polished consensus sequence file and all intermediate files will be written to the designated output path. User must specify threads to use, and a medaka model for the basecaller used.

# Example usage: ./medaka_polish.sh -r read.fastq -d draft.fa -o outdir/ -t 2 -m r941_min_sup_g507

set -eo pipefail

### FUNCTIONS ###

# prints help function
function usage {
	echo -e "\nThis is a wrapper script that automates four rounds of medaka polishing workflow." 
	echo -e "\nFlags:" 
	echo -e "-r \t Input fastq file containing the nanopore reads."
	echo -e "-d \t Input fasta file containing the draft sequence to be polished. Required name format: <sample_name>.draft.fa"
	echo -e "-o \t Output path."
	echo -e "-t \t Number of threads."
	echo -e "-m \t Medaka model to use. See Medaka documentation for details."
	echo -e "\nExample usage: \n ./medaka_polish.sh -r reads.fastq -d draft.fa -o results/ -t 2 -m r941_min_sup_g507"
	exit 1
}

# parses command line flags
function getoptsfunction {
	local OPTIND
	while getopts ":r:d:o:t:m:h" flag; do
		case $flag in
			r) READS=$OPTARG;;
			d) DRAFT=$OPTARG 
				echo ""; echo "Draft sequence used: ${DRAFT}" ;;
			o) OUTPUT=$OPTARG ;;
			t) THREADS=$OPTARG ;;
			m) MODEL=$OPTARG ;;
			h) usage ;;
			:) echo "Argument missing from -${OPTARG} option" >&2; exit 1 ;;
			\?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
		esac
	done
	shift $(( OPTIND-1 )) # Tells getopts to move onto the next argument
	
	if [[ -z "$READS" ]] ; then echo "Missing input read file (-r)."; usage; >&2; exit 1; fi
	if [[ -z "$DRAFT" ]] ; then echo "Missing input draft sequence file (-d)."; usage; >&2; exit 1; fi
	if [[ -z "$OUTPUT" ]] ; then echo "Missing output dir path (-o)."; usage; >&2; exit 1; fi
	if [[ -z "$THREADS" ]] ; then echo "Missing number of threads to use (-t)."; usage; >&2; exit 1; fi
	if [[ -z "$MODEL" ]] ; then echo "Missing medaka model to use (-m)."; usage; >&2; exit 1; fi

}

getoptsfunction "$@"

# prints polishing round report
function polishing_round_report() {
	echo ""
	echo "Polishing round $1/4..."
	echo ""
} 


### CODES ###

mkdir -p ${OUTPUT}
echo "" ; echo "###### Starting medaka polishing ######"

# extract sample name from draft filename
SAMPLE=$( basename ${DRAFT} )
SAMPLE=${SAMPLE%.draft.fa}

# run four polishing rounds
for i in $( seq 1 4 ); do

	polishing_round_report $i
	
	# filenames handling
	if [[ $i = 1 ]] # round 1 uses draft
	then
		mini_align_input=${DRAFT}
		prev_prefix=${OUTPUT}/${SAMPLE}.draft
	else # rounds starting from 2 use previous round
		prev=$( expr $i - 1 )
		mini_align_input=${OUTPUT}/${SAMPLE}.polish.round${prev}.fasta
		prev_prefix=${OUTPUT}/${SAMPLE}.polish.round${prev}
	fi
	PREFIX=${OUTPUT}/${SAMPLE}.polish.round$i
	
	# run polishing pipeline
	mini_align -i ${READS} -r ${mini_align_input} -m -p ${prev_prefix} -t ${THREADS}
	medaka consensus ${prev_prefix}.bam ${PREFIX}.hdf --model ${MODEL} --chunk_len 6000 --threads ${THREADS}
	medaka stitch ${PREFIX}.hdf ${mini_align_input} ${PREFIX}.fasta

done


# produce final consensus and store intermediate files
POLISH_ROUND4_FASTA=${OUTPUT}/${SAMPLE}.polish.round4.fasta
FINAL_CONSENSUS_FASTA=${OUTPUT}/${SAMPLE}.final_consensus.fa


if [[ -f $POLISH_ROUND4_FASTA ]]

then
	cp $POLISH_ROUND4_FASTA $FINAL_CONSENSUS_FASTA
	
	# change header from read id (centroid) to sample name
	sed -i "1s/^.*/>${SAMPLE}/" $FINAL_CONSENSUS_FASTA
	
	mkdir -p ${OUTPUT}/intermediate_files
	for file in ${OUTPUT}/*; do
		if ! { [[ "$file" == *intermediate_files ]] || [[ "$file" == *final* ]]; }
		then
			mv $file ${OUTPUT}/intermediate_files/.
		fi
	done
	
	echo ""
	echo "The final polished consensus sequence is written to: $FINAL_CONSENSUS_FASTA"
	echo ""
	echo "Medaka polishing done!"

else
	echo "Polishing round4 output does not exist. Please check ${SAMPLE}.medaka_polish.log in logs directory for errors."
	
fi
