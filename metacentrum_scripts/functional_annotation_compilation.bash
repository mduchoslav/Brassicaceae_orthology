### Metacentrum script

#PBS -N functional_annotation_compilation
#PBS -l select=1:ncpus=1:mem=32gb:scratch_local=200gb
#PBS -l walltime=24:00:00 
#PBS -m ae

## It is needed to set the $species variable during submitting the job like this:
# for species in Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Cardamine_glauca Noccaea_praecox
# do
# echo "Submitting job for species: $species"
# qsub  -v "species=$species" InterProScan_all_species.bash
# done

echo "Species for this job: $species" | ts '[%Y-%m-%d %H:%M:%S]'

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# move to working dir
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3

# load R
module load r/4.1.3-gcc-10.2.1-6xt26dl

# run R script
Rscript --verbose "functional_annotation_compilation_executable.r" $species || { echo >&2 "Rscript failed!"; exit 1; }

echo "Run ended: $species" | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch