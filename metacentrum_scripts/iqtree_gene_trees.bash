### Metacentrum script

#!/bin/bash
#PBS -N iqtree_gene_trees
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=10gb
#PBS -l walltime=24:00:00
#PBS -m ae

set -euo pipefail

############################
# User-defined variables
############################

# Directory with MAFFT alignments
INPUT_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/mafft_alignments"

# File containing selected orthogroup IDs, one per line
ORTHOGROUP_LIST="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/input_orthogroups/complete_orthogroups_max30.ids.txt"

# Directory for IQ-TREE results
OUTPUT_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/iqtree_gene_trees"

# Number of parallel IQ-TREE jobs
NCPU=16

# Batch size
BATCH_SIZE=500

# Batch number.
# Can be supplied using:
# qsub -v BATCH=3 iqtree_gene_trees.sh
#
# If BATCH is not supplied, batch 1 is used.
BATCH="${BATCH:-1}"

############################
# Calculate batch range
############################

START=$(( (BATCH - 1) * BATCH_SIZE + 1 ))
END=$(( BATCH * BATCH_SIZE ))

############################
# Environment
############################

module load iqtree

mkdir -p "$OUTPUT_DIR"

cd "$OUTPUT_DIR"

echo "Starting IQ-TREE gene tree inference"
echo "Input directory:    $INPUT_DIR"
echo "Orthogroup list:    $ORTHOGROUP_LIST"
echo "Output directory:   $OUTPUT_DIR"
echo "Batch:              $BATCH"
echo "Batch range:        $START-$END"
echo "Parallel jobs:      $NCPU"
echo "Start time:         $(date "+%Y-%m-%d %H:%M:%S")"
echo "IQ-TREE version:"
iqtree3-mpi --version 2>&1
echo

############################
# Prepare batch list
############################

BATCH_LIST="${SCRATCHDIR}/orthogroups_batch_${BATCH}.txt"

grep -v '^[[:space:]]*$' "$ORTHOGROUP_LIST" |
    sed -n "${START},${END}p" \
    > "$BATCH_LIST"

N_SELECTED=$(wc -l < "$BATCH_LIST")

echo "Orthogroups in this batch: $N_SELECTED"
echo

if [[ "$N_SELECTED" -eq 0 ]]; then
    echo "No orthogroups found for batch $BATCH."
    clean_scratch
    exit 0
fi

############################
# Gene tree function
############################

tree_one() {

    og="$1"

    # Remove possible Windows carriage returns
    og="${og//$'\r'/}"

    # Skip empty lines
    [[ -z "$og" ]] && return

    infile="${INPUT_DIR}/${og}.aln.fa"
    prefix="${OUTPUT_DIR}/${og}"

    # Check whether the input alignment exists
    if [[ ! -f "$infile" ]]; then
        echo "WARNING: Input alignment not found: $infile" >&2
        return
    fi

    # Skip already completed trees
    if [[ -s "${prefix}.treefile" ]]; then
        echo "Skipping existing tree: $og"
        return
    fi

    echo "Inferring tree: $og"

    iqtree3-mpi \
        -s "$infile" \
        -m LG+G4 \
        -T 1 \
        --prefix "$prefix" \
        --quiet
}

export -f tree_one
export INPUT_DIR
export OUTPUT_DIR

############################
# Run orthogroups in batch
############################

xargs -n 1 -P "$NCPU" \
    bash -c 'tree_one "$1"' _ \
    < "$BATCH_LIST"

############################
# Summary
############################

N_OUTPUT=0

while read -r og; do

    og="${og//$'\r'/}"

    if [[ -s "${OUTPUT_DIR}/${og}.treefile" ]]; then
        ((N_OUTPUT+=1))
    fi

done < "$BATCH_LIST"

echo
echo "IQ-TREE gene tree inference finished"
echo "Batch:                 $BATCH"
echo "Orthogroups in batch:  $N_SELECTED"
echo "Completed trees:       $N_OUTPUT"
echo "End time:              $(date "+%Y-%m-%d %H:%M:%S")"

clean_scratch