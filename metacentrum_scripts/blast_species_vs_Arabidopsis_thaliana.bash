### Metacentrum script

#PBS -N blast_species_vs_Arabidopsis_thaliana
#PBS -l select=1:ncpus=8:mem=8gb:scratch_local=20gb
#PBS -l walltime=2:00:00 
#PBS -m ae

## It is needed to set the $species variable during submitting the job like this:
# for species in Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Cardamine_glauca Noccaea_praecox
# do
# echo "Submitting job for species: $species"
# qsub  -v "species=$species" blast_species_vs_Arabidopsis_thaliana.bash
# done

# define variables
query_fasta="/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/primary_transcripts/$species.fasta"
target_fasta=/storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/primary_transcripts/Arabidopsis_thaliana.fasta
# outfmt=0
outfmt='7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen'
# -outfmt 0 = Pairwise
# -outfmt 6 = Tabular
# -outfmt 7 = Tabular with comment lines
output_file=$species"_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv"
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

# Make blast database
makeblastdb -in $target_fasta_short -dbtype prot -out target_db

echo "Database done."  | ts '[%Y-%m-%d %H:%M:%S]'

# run blast search
blastp -evalue 0.1 -db target_db -query $query_fasta_short -out $output_file -outfmt "$outfmt" -num_threads 8

echo "Blast search done." | ts '[%Y-%m-%d %H:%M:%S]'

# move the output to user's DATADIR or exit in case of failure
mkdir -p $output_dir # will make folder if it doesn't exist
cp $output_file $output_dir/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Copying output file $output_file to $output_dir done." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch