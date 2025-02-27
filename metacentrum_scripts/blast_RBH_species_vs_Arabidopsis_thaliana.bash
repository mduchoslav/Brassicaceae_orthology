### Metacentrum script

#PBS -N blast_RBH_species_vs_Arabidopsis_thaliana
#PBS -l select=1:ncpus=8:mem=4gb:scratch_local=20gb
#PBS -l walltime=2:00:00 
#PBS -m ae

## It is needed to set the $species variable during submitting the job like this:
# for species in Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Cardamine_glauca Noccaea_praecox
# do
# echo "Submitting job for species: $species"
# qsub  -v "species=$species" blast_RBH_species_vs_Arabidopsis_thaliana.bash
# done

# define variables
# It is reciprocal here, so it does not matter which is query and which target.
query_fasta="/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/primary_transcripts/$species.fasta"
target_fasta=/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/primary_transcripts/Arabidopsis_thaliana.fasta
output_file=$species"_vs_Arabidopsis_thaliana_rbh-i70-c50.txt"
output_dir=/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/blast_results


# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy input files to scratch directory
# if the copy operation fails, issue error message and exit
cp -t $SCRATCHDIR $query_fasta $target_fasta || { echo >&2 "Error while copying input file(s)!"; exit 2; }

echo "Input files copied."  | ts '[%Y-%m-%d %H:%M:%S]'

# move into scratch directory
cd $SCRATCHDIR 

# fasta files without pathes
query_fasta_short=$(echo $query_fasta | sed "s,/.*/,,g")
target_fasta_short=$(echo $target_fasta | sed "s,/.*/,,g")

# load module
module load blast-plus/2.16.0-gcc-10.2.1-bgzrrrz
module load python36-modules # needed just when --nr option is used

# run the Reciprocal blast hit search
python3 /storage/brno12-cerit/home/duchmil/SW/blast_rbh/blast_rbh.py -a prot -t blastp -o $output_file --threads=8 --nr $query_fasta_short $target_fasta_short
# also --nr parameter in one variant (filter to only one copy of identical sequences)

echo "blast_rbh.py script finished." | ts '[%Y-%m-%d %H:%M:%S]'

# move the output to user's DATADIR or exit in case of failure
mkdir -p $output_dir # will make folder if it doesn't exist
cp $output_file $output_dir/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Copying output file $output_file to $output_dir done." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch