### Metacentrum script

#PBS -N InterProScan_all_species
#PBS -l select=1:ncpus=8:mem=96gb:scratch_local=200gb
#PBS -l walltime=24:00:00 
#PBS -m ae

## It is needed to set the $species variable during submitting the job like this:
# for species in Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Cardamine_glauca Noccaea_praecox
# do
# echo "Submitting job for species: $species"
# qsub  -v "species=$species" InterProScan_all_species.bash
# done

echo "Species for this job: $species"

# define variables
query_fasta="/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/interproscan/primary_transcripts_without_stops/$species.fasta"
output_dir=/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/interproscan

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# move into scratch directory
cd $SCRATCHDIR 

# load Java
module load openjdk

# run InterProScan
/storage/brno12-cerit/home/duchmil/SW/InterProScan/interproscan-5.72-103.0/interproscan.sh --goterms --pathways --enable-tsv-residue-annot --tempdir $SCRATCHDIR --cpu 8 -i $query_fasta --output-file-base $species

echo "Main calculation done." | ts '[%Y-%m-%d %H:%M:%S]'

# compress output files
pigz -v -p 8 $species.*

# move the output to user's DATADIR or exit in case of failure
mkdir -p $output_dir # will make folder if it doesn't exist
cp -v $species.* $output_dir/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Copying output file to $output_dir done." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch