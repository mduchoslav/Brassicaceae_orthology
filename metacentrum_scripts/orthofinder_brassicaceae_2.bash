### Metacentrum script

#PBS -N orthofinder_brassicaceae_2
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=200gb
#PBS -l walltime=23:50:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

# loading Python and the modules needed
module add python36-modules-gcc

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy Orthofinder
cp -r /storage/brno12-cerit/home/duchmil/orthofinder/OrthoFinder_source $SCRATCHDIR || { echo >&2 "Error while copying Orthofinder file(s)!"; exit 2; }

# copy input files to scratch directory
# if the copy operation fails, issue error message and exit
cp -r $DATADIR/primary_transcripts $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

# move into scratch directory
cd $SCRATCHDIR 

# running OrthoFinder
$SCRATCHDIR/OrthoFinder_source/orthofinder.py -y -t 8 -f $SCRATCHDIR/primary_transcripts/ -o orthofinder_results -n brassicaceae_2

# move the output to user's DATADIR or exit in case of failure
cp -r orthofinder_results $DATADIR/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch