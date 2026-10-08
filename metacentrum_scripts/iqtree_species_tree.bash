### Metacentrum script

#!/bin/bash
#PBS -N iqtree_species_tree
#PBS -l select=1:ncpus=32:mem=128gb:scratch_local=100gb
#PBS -l walltime=168:00:00
#PBS -m ae

set -euo pipefail

############################
# User-defined variables
############################

# Directory with trimmed locus alignments
ALIGN_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/trimmed_alignments_2"

# Output directory
OUTPUT_DIR="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/phylo_tree/iqtree_species_tree"

# Output prefix
PREFIX="brassicaceae_species_tree"

# Number of CPUs
NCPU=32

############################
# Environment
############################

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# move into scratch directory
cd $SCRATCHDIR

echo "Start copying input data." | ts '[%Y-%m-%d %H:%M:%S]'

# copy input data
mkdir iqtree_input
cp ${ALIGN_DIR}/* iqtree_input

echo "Input data copying done." | ts '[%Y-%m-%d %H:%M:%S]'

# dir for output
mkdir iqtree_output

# load SW
module load iqtree

echo "Starting IQ-TREE species tree inference"
echo "Alignment directory: $ALIGN_DIR"
echo "Output directory:    $OUTPUT_DIR"
echo "CPUs:                $NCPU"
echo "Start time:          $(date "+%Y-%m-%d %H:%M:%S")"
echo
echo "IQ-TREE version:"
iqtree3-mpi --version 2>&1
echo

############################
# Run IQ-TREE
############################

iqtree3-mpi \
    -p "${SCRATCHDIR}/iqtree_input" \
    --seqtype AA \
    -m MFP+MERGE \
    -rcluster 10 \
    -B 1000 \
    -alrt 1000 \
    -T "$NCPU" \
    --prefix "${SCRATCHDIR}/iqtree_output/${PREFIX}"

echo "Start copying output data." | ts '[%Y-%m-%d %H:%M:%S]'

# copy output data
mkdir -p "$OUTPUT_DIR"
cp -R ${SCRATCHDIR}/iqtree_output/* $OUTPUT_DIR

echo "Output data copying done." | ts '[%Y-%m-%d %H:%M:%S]'

############################
# Summary
############################

echo
echo "IQ-TREE analysis finished"
echo "End time: $(date "+%Y-%m-%d %H:%M:%S")"
echo
echo "Main output tree:"
echo "  ${OUTPUT_DIR}/${PREFIX}.treefile"

clean_scratch