### Metacentrum script

#!/bin/bash
#PBS -N mafft_orthogroups
#PBS -l select=1:ncpus=16:mem=2gb:scratch_local=10gb
#PBS -l walltime=1:00:00
#PBS -m ae

set -euo pipefail

############################
# User-defined variables
############################

# Directory with OrthoFinder orthogroup FASTA files
INPUT_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/input_orthogroups/random_single_copy_fastas"

# File containing selected orthogroup IDs, one per line
ORTHOGROUP_LIST="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/input_orthogroups/complete_orthogroups_max30.ids.txt"

# Directory for MAFFT alignments
OUTPUT_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/mafft_alignments"

# Number of parallel MAFFT jobs
NCPU=16

############################
# Environment
############################

module load mafft

mkdir -p "$OUTPUT_DIR"

cd "$OUTPUT_DIR"

echo "Starting MAFFT alignments"
echo "Input directory:    $INPUT_DIR"
echo "Orthogroup list:    $ORTHOGROUP_LIST"
echo "Output directory:   $OUTPUT_DIR"
echo "Parallel jobs:      $NCPU"
echo "Start time:         $(date "+%Y-%m-%d %H:%M:%S")"
echo "MAFFT version:"
mafft --version 2>&1
echo

############################
# Alignment function
############################

align_one() {

    og="$1"

    # Remove possible Windows carriage returns
    og="${og//$'\r'/}"

    # Skip empty lines
    [[ -z "$og" ]] && return

    infile="${INPUT_DIR}/${og}.fa"
    outfile="${OUTPUT_DIR}/${og}.aln.fa"

    # Check whether the input file exists
    if [[ ! -f "$infile" ]]; then
        echo "WARNING: Input file not found: $infile" >&2
        return
    fi

    # Skip already completed alignments
    if [[ -s "$outfile" ]]; then
        echo "Skipping existing alignment: $og"
        return
    fi

    echo "Aligning: $og"

    mafft \
        --auto \
        --thread 1 \
        "$infile" \
        > "$outfile"
}

export -f align_one
export INPUT_DIR
export OUTPUT_DIR

############################
# Run selected orthogroups
############################

grep -v '^[[:space:]]*$' "$ORTHOGROUP_LIST" |
    xargs -n 1 -P "$NCPU" bash -c 'align_one "$1"' _

############################
# Summary
############################

N_SELECTED=$(grep -cv '^[[:space:]]*$' "$ORTHOGROUP_LIST")

N_OUTPUT=$(find "$OUTPUT_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.aln.fa" \
    -size +0c |
    wc -l)

echo
echo "MAFFT alignment finished"
echo "Selected orthogroups: $N_SELECTED"
echo "Completed alignments: $N_OUTPUT"
echo "End time:             $(date "+%Y-%m-%d %H:%M:%S")"

# clean the SCRATCH directory
clean_scratch