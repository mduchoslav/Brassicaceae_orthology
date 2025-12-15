Enriching list of orthologues
================
Miloš Duchoslav
2025-11

- [Introduction](#introduction)
- [Strategy](#strategy)
- [Length of proteins](#length-of-proteins)
  - [Calculation of protein length from
    fasta](#calculation-of-protein-length-from-fasta)
- [Running BLAST and reciprocal
  BLAST](#running-blast-and-reciprocal-blast)
  - [BLAST](#blast)
  - [Reciprocal best hit (RBH) BLAST](#reciprocal-best-hit-rbh-blast)
- [Enriching list of orthologues - compilation of
  results](#enriching-list-of-orthologues---compilation-of-results)
  - [Load functions etc.](#load-functions-etc)
  - [Reading common files](#reading-common-files)
  - [Big loop](#big-loop)
- [Explanation of columns in the output
  table](#explanation-of-columns-in-the-output-table)
- [Improving table with stats](#improving-table-with-stats)
- [End](#end)

# Introduction

This RMarkdown file (or its markdown version for GitHub) documents
supplementing the list of orthologues from OrthoFinder with results from
BLAST and with orthogroups from OrthoFinder to get *Arabidopsis
thaliana* homologues for as many genes of the target species as
possible. The homologues will be used in other scripts to transfer
functional annotation from *Arabidopsis thaliana* genes.

This file includes:

1.  BASH code that I ran at MetaCentrum (Czech national grid
    infrastructure) running PBS scheduling system for batch jobs.
2.  R code that I ran locally.

### SW installation and versions

The SW installation instructions and versions of SW used is described in
[Installation_of_SW.md](Installation_of_SW.md).

# Strategy

### Sources of information

1.  Orthologues tables from OrthoFinder
    - Tables from `Orthologues` folder from OrthoFinder results.
2.  Phylogenetic hierarchical orthogroups from OrthoFinder (N0 table)
    - Table `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` from
      OrthoFinder results.
    - It is more finely resolved compared to Orthogroups table.
3.  Orthogroups table from OrthoFinder
    - Table `Orthogroups/Orthogroups.tsv` from OrthoFinder results.
    - The orthogroups in this table are broader compared to N0
      orthogroups.
4.  Reciprocal best hit (RBH) BLAST
5.  Simple BLAST

The sources other than Orthologues tables from OrthoFinder might not
give us the real orthologues in *Arabidopsis thaliana*, but possibly
other kinds of homologues. However, these homologues will likely have
similar functions, so it will be useful for transferring of functional
annotation from these genes.

### Final selection of orthologues

Strategy for column with final orthologues (column
`Arabidopsis_thaliana`):

1.  Use orthologues from orthologues table.
2.  If missing, use Blast RBH.
3.  If missing, use genes from N0 hierarchical orthogroups.
4.  If missing, use genes from broad orthogroups.
5.  If missing, use genes hits from simple Blast (but only those with
    BLAST_pident \> 40 & BLAST_qcovhsp \> 50)
    - BLAST_pident: Percentage of identical matches (in local alignment)
    - BLAST_qcovhsp: Query coverage per HSP (%)
      - HSP = High-scoring Segment Pair (local alignment with no gaps)

Strategy for column with single orthologues (column
`single_Arabidopsis_thaliana`):

1.  If there is only one *A. thaliana* gene in column
    `Arabidopsis_thaliana`, use that one.
2.  Else, if there is one and only RBH and it is among genes in column
    `Arabidopsis_thaliana`, take the RBH.
3.  Else, take the gene from the column `Arabidopsis_thaliana` that had
    the highest bitscore in the simple BLAST results. If several have
    the highest bitscore, take the first.
4.  If the genes in column `Arabidopsis_thaliana` are not among BLAST
    hits, take nothing.

**Recommendation**: It is appealing to use the “single” orthologues
(`single_Arabidopsis_thaliana`), because it makes things easier.
However, I recommend to do that only in cases when it is absolutely
necessary and otherwise use the column where there are multiple
orthologues in some cases (`Arabidopsis_thaliana`). The real orthology
is not always one-to-one due to multiplications of genes in some species
(there is a good explanation in [OrthoFinder
GitHub](https://github.com/davidemms/OrthoFinder)). The column
`Arabidopsis_thaliana` should better reflect the real biology.

# Length of proteins

I will add length of protein sequence (of the longest transcript) to the
list of genes. It is good for comparison with length of obtained
homologues. Also, it will help me to get all genes in the list.

## Calculation of protein length from fasta

``` sh
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/
mkdir -p protein_length_primary_transcripts

for species_file in primary_transcripts/*.fasta
do
echo "Calculating length for: $species_file"
species=$(sed -E 's-(^.*/)|(\.fasta)--g' <(echo $species_file))
# echo $species

## Calculate length of fasta sequences
# Script from AI, tested
awk '
BEGIN { name=""; seq="" }
/^>/ {
    if (name != "") {
        seq_clean = seq
        sub(/[ *]+$/, "", seq_clean)   # remove trailing * only
        print name "\t" length(seq_clean)
    }
    # extract first word of header without ">"
    header = substr($0, 2)
    split(header, a, /[ \t]/)
    name = a[1]
    seq = ""
    next
}
{
    line = $0
    sub(/[ \t\r\n]+$/, "", line)
    seq = seq line
}
END {
    if (name != "") {
        seq_clean = seq
        sub(/[ *]+$/, "", seq_clean)
        print name "\t" length(seq_clean)
    }
}
' $species_file > protein_length_primary_transcripts/$species"_protein_length.tsv"

done
```

# Running BLAST and reciprocal BLAST

Approach similar to:  
\> Bray, Sian M., Tuomas Hämälä, Min Zhou, Silvia Busoms, Sina Fischer,
Stuart D. Desjardins, Terezie Mandáková, et al. “Kinetochore and Ionomic
Adaptation to Whole-Genome Duplication in Cochlearia Shows Evolutionary
Convergence in Three Autopolyploids.” Cell Reports 43, no. 8 (August 27,
2024). <https://doi.org/10.1016/j.celrep.2024.114576>.

The code for the paper is at Sian Bray’s
[GitHub](https://github.com/Sian-Bray/Cochlearia_2024).

## BLAST

### Script for running Blast on Metacentrum

``` bash
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
query_fasta="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/primary_transcripts/$species.fasta"
target_fasta="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/primary_transcripts/Arabidopsis_thaliana.fasta"
# outfmt=0
outfmt='7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen'
# -outfmt 0 = Pairwise
# -outfmt 6 = Tabular
# -outfmt 7 = Tabular with comment lines
output_file=$species"_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv"
output_dir=/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/blast_results


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

# Resources used: Up to 46 min (mostly 11-20 min), 64-99% CPU, up to 1 GB memory.
```

### Submitting jobs for several species

``` bash
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/metacentrum_scripts

for species in Aethionema_saxatile Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Arabidopsis_lyrata_Rawat Arabis_alpina Brassica_napus Brassica_oleracea Brassica_rapa Camelina_sativa Capsella_rubella Cardamine_amara Cardamine_glauca Cardamine_hirsuta Cochlearia_excelsa Conringia_planisiliqua Erysimum_linariifolium Euclidium_syriacum Eutrema_salsugineum Noccaea_praecox Odontarrhena_muralis Raphanus_sativus
do
echo "Submitting job for species: $species"
qsub  -v "species=$species" blast_species_vs_Arabidopsis_thaliana.bash
done
```

### Compressing results

``` bash
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/blast_results
# compress the result
gzip *_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv
```

### Format of the output file

`outfmt='7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen'`

| Column | NCBI name | Description                                  |
|--------|-----------|----------------------------------------------|
| 1      | qaccver   | Query accession and version                  |
| 2      | saccver   | Subject accession and version                |
| 3      | pident    | Percentage of identical matches              |
| 4      | length    | Alignment length                             |
| 5      | mismatch  | Number of mismatches                         |
| 6      | gapopen   | Number of gap openings                       |
| 7      | qstart    | Start of alignment in query                  |
| 8      | qend      | End of alignment in query                    |
| 9      | sstart    | Start of alignment in subject (database hit) |
| 10     | send      | End of alignment in subject (database hit)   |
| 11     | evalue    | Expectation value (E-value)                  |
| 12     | bitscore  | Bit score                                    |
| 13     | qcovhsp   | Query coverage per HSP (%)                   |
| 14     | qlen      | Query sequence length                        |
| 15     | slen      | Subject sequence length                      |

Suggestions from
<https://github.com/peterjc/galaxy_blast/blob/master/tools/blast_rbh/blast_rbh.xml>:

    If you are trying to use BLAST RBH matches to identify candidate orthologues
    or transfer annotation, you *must* use a percentage identity and minimum
    coverage threshold or similiar. See:

    Punta and Ofran (2008) The Rough Guide to In Silico Function Prediction,
    or How To Use Sequence and Structure Information To Predict Protein
    Function. PLoS Comput Biol 4(10): e1000160.
    https://doi.org/10.1371/journal.pcbi.1000160

    The defaults are to require 70% sequence identity over the aligned region
    (using ``pident`` in the BLAST+ tabular output), and that the HSP alignment
    covers at least 50% of the query sequence (using ``qcovhsp`` in the BLAST+
    tabular output).

## Reciprocal best hit (RBH) BLAST

### Running RBH on Metacentrum

``` bash
### Metacentrum script

#PBS -N blast_RBH_species_vs_Arabidopsis_thaliana
#PBS -l select=1:ncpus=8:mem=4gb:scratch_local=20gb
#PBS -l walltime=4:00:00 
#PBS -m ae

## It is needed to set the $species variable during submitting the job like this:
# for species in Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Cardamine_glauca Noccaea_praecox
# do
# echo "Submitting job for species: $species"
# qsub  -v "species=$species" blast_RBH_species_vs_Arabidopsis_thaliana.bash
# done

# define variables
# It is reciprocal here, so it does not matter which is query and which target.
query_fasta="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/primary_transcripts/$species.fasta"
target_fasta="/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/primary_transcripts/Arabidopsis_thaliana.fasta"
output_file=$species"_vs_Arabidopsis_thaliana_rbh-i70-c50.txt"
output_dir=/storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/blast_results


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

# Resources: up to 1 h, 68-99 % CPU, up to 660 MB memory
```

### Submitting jobs for several species

``` bash
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/metacentrum_scripts

for species in Aethionema_saxatile Alyssum_gmelinii Arabidopsis_arenosa Arabidopsis_lyrata_NCBI Arabidopsis_lyrata_Rawat Arabis_alpina Brassica_napus Brassica_oleracea Brassica_rapa Camelina_sativa Capsella_rubella Cardamine_amara Cardamine_glauca Cardamine_hirsuta Cochlearia_excelsa Conringia_planisiliqua Erysimum_linariifolium Euclidium_syriacum Eutrema_salsugineum Noccaea_praecox Odontarrhena_muralis Raphanus_sativus
do
echo "Submitting job for species: $species"
qsub  -v "species=$species" blast_RBH_species_vs_Arabidopsis_thaliana.bash
done
```

``` bash
cd /storage/brno12-cerit/home/duchmil/Brassicaceae_orthology/brassicaceae_3/blast_results
wc -l *_vs_Arabidopsis_thaliana_rbh-i70-c50.txt


grep ';' Arabidopsis_arenosa_vs_Arabidopsis_thaliana_rbh-i70-c50.txt
# Genes with identical protein sequences are merged, the gene IDs are separated by ';'.
```

Genes with identical protein sequences are merged, the gene IDs are
separated by ‘;’.

Number of RBHs: 15596
Alyssum_gmelinii_vs_Arabidopsis_thaliana_rbh-i70-c50.txt 21192
Arabidopsis_arenosa_vs_Arabidopsis_thaliana_rbh-i70-c50.txt 22534
Arabidopsis_lyrata_NCBI_vs_Arabidopsis_thaliana_rbh-i70-c50.txt 18265
Cardamine_glauca_vs_Arabidopsis_thaliana_rbh-i70-c50.txt 17388
Noccaea_praecox_vs_Arabidopsis_thaliana_rbh-i70-c50.txt

There was still this error message:
`Warning: Sequences with tied best hits found, you may have duplicates/clusters`
I guess it might mean that there are some proteins that have identical
sequence for some domain and the hit was to that domain.

# Enriching list of orthologues - compilation of results

R script for compilation of results from the different sources.

## Load functions etc.

``` r
setwd("D:/!ecolgen/Brassicaceae_orthology/brassicaceae_3")
old.par<-par(no.readonly = T)

### Function to split multiple values in one row to several rows 
### - fast version based on data.table package

# Just one column will be splitted, the other will be repeated.
# x - data.frame to be splitted
# col.split - number or name of the column where the values should be splitted
# split - separator of the values to be splitted

library(data.table)

df_split_dt <- function(x, col.split, split = ", ") {
  dt <- as.data.table(x)
  if(is.numeric(col.split)) {
    colname <- names(dt)[col.split]
  } else {
    colname <- col.split
  }
  dt[, strsplit(get(colname), split), by = setdiff(names(dt), colname)][, setnames(.SD, "V1", colname) ][, names(dt), with = F]
}


### Function for printing messages with time stamps

# Arguments:
# ... - Anything that will be pasted after the time stamp.
# sep - Separator of the arguments printed after the time stamp. Defaults to space.

t.print <- function(..., sep = " ") {
    print(paste0("[", format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "] ", paste(..., sep = sep)))
  }
```

## Reading common files

``` r
o.groups_a <- read.table(file = "orthofinder_results/Results_brassicaceae_3/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", sep = "\t", header = T)

colnames(o.groups_a)

o.groups_b <- read.table(file = "orthofinder_results/Results_brassicaceae_3/Orthogroups/Orthogroups.tsv", sep = "\t", header = T)
```

## Big loop

This loop cycles through species where the orthologues from Orthofinder
should be supplemented. Besides supplementing orthologues it collects
statistics.

``` r
# read species names
blast.files <- list.files(path = "blast_results/", pattern = ".tsv.gz")
species <- sub(pattern = "_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv.gz", replacement = "", x = blast.files, fixed = T)

# chosen species for testing of the loop
one.species <- "Aethionema_saxatile"
one.species <- "Alyssum_gmelinii"
one.species <- "Arabidopsis_arenosa"

# Prepare dataframe for statistics
stats.ortho <- data.frame()


### Start of the big loop

for(one.species in species) {
  
  # Save console output to log file
    dir.create(path = "supplemented_orthologues/logs", recursive = T, showWarnings = F)
    log.file <- file(paste0("supplemented_orthologues/logs/", one.species, "_supplementing_orthologues.log"), open = "wt")
  sink(file = log.file, split = T)
  sink(file = log.file, type = "message")
  
  ### Lenght of proteins
  # This is also to get IDs of all protein-coding genes to start with
  prot.length <- fread(file = paste0("protein_length_primary_transcripts/", one.species, "_protein_length.tsv"), 
                                         header = F, sep = "\t", na.strings = "",
                                         # comment.char = "", 
                                         quote = "")
  colnames(prot.length) <- c(one.species, "prot_length")

  
  ### Orthologues table from Orthofinder
  
  t.print("Reading orthologues table for", one.species)
  
  # ortho1 <- read.table(file = paste0("orthofinder_results/Results_brassicaceae_3/Orthologues/Orthologues_", 
  #                                                             one.species, "/", one.species, "__v__Arabidopsis_thaliana.tsv"), 
  #                                 sep = "\t", header = T)
  ortho1 <- fread(file = paste0("orthofinder_results/Results_brassicaceae_3/Orthologues/Orthologues_", 
                                     one.species, "/", one.species, "__v__Arabidopsis_thaliana.tsv"), 
                       sep = "\t", header = T)
  # remove "Orthogroup" column, it is same as in `Orthogroups/Orthogroups.tsv`
  # ortho1b <- ortho1[, c(one.species, "Arabidopsis_thaliana")]
  ortho1b <- ortho1[, c(one.species, "Arabidopsis_thaliana"), with = F]
  
  # rename the column to show that it comes from orthologues (OL) table
  colnames(ortho1b)[2] <- paste0("OL_", colnames(ortho1b)[2])
  
  # splitting lines to have one gene of the desired species per line
  ortho.split <- df_split_dt(x = ortho1b, col.split = one.species, split = ", ")
  
  ## putting in big table
  ortho.big.0 <- merge(x = prot.length, y = ortho.split, by = one.species, all = T)

  
  
  ### Supplementing orthologues from `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`
  
  t.print("Supplementing orthologues from `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` for", one.species)
  
  # choosing right columns
  # (the df.split function is then much faster)
  # Also OG is not needed, it is always the same as in `Orthogroups/Orthogroups.tsv`
  o.groups_a_2 <- o.groups_a[, c("HOG", one.species, "Arabidopsis_thaliana")]
  colnames(o.groups_a_2) <- c("N0_HOG", one.species, "N0_Arabidopsis_thaliana")
  
  # rows with both target species and thaliana genes
  o.groups.a.t_a <- o.groups_a_2[o.groups_a_2[, one.species] != "" & o.groups_a_2$N0_Arabidopsis_thaliana != "", ]
  
  # splitting lines to have one gene of the desired species per line
  o.groups.split_a <- df_split_dt(x = o.groups.a.t_a, col.split = one.species, split = ", ")
  
  ## putting in big table
  ortho.big.1 <- merge(x = ortho.big.0, y = o.groups.split_a, by = one.species, all = T)
  
  
  
  ### Supplementing orthologues from `Orthogroups/Orthogroups.tsv`
  
  t.print("Supplementing orthologues from `Orthogroups/Orthogroups.tsv` for", one.species)
  
  # choosing right columns
  # (the df.split function is then much faster)
  o.groups_b_2 <- o.groups_b[, c("Orthogroup", one.species, "Arabidopsis_thaliana")]
  colnames(o.groups_b_2) <- c("OG_Orthogroup", one.species, "OG_Arabidopsis_thaliana")
  # colnames(o.groups_b)
  # colnames(o.groups_b_2)
  
  # rows with both target species and thaliana genes
  o.groups.a.t_b <- o.groups_b_2[o.groups_b_2[, one.species] != "" & o.groups_b_2$OG_Arabidopsis_thaliana != "", ]
  
  # splitting lines to have one gene of the desired species per line
  o.groups.split_b <- df_split_dt(x = o.groups.a.t_b, col.split = which(colnames(o.groups.a.t_b) == one.species), split = ", ")
  
  
  ## putting in big table
  ortho.big.2 <- merge(x = ortho.big.1, y = o.groups.split_b, by = one.species, all = T)
  
  
  ### Merging with blast RBH table
  
  t.print("Merging with blast RBH table for", one.species)
  
  # read Blast reciprocal best hits (RBH)
  rbh.pre <- fread(file = paste0("blast_results/", one.species, "_vs_Arabidopsis_thaliana_rbh-i70-c50.txt"), sep = "\t", 
                                 header = T, 
                                 # comment.char = ""
                                 )
  
  colnames(rbh.pre)[1] <- one.species
  colnames(rbh.pre)[-1] <- paste0("RBH_", colnames(rbh.pre)[-1])
  
  # splitting RBH
  rbh <- df_split_dt(x = rbh.pre, col.split = 1, split = ";")
  
  # change the separator for thaliana genes in RBH table (to be the same as in other tables)
  rbh$RBH_B_id <- gsub(pattern = ";", replacement = ", ", x = rbh$RBH_B_id)
  
  
  # merging to large table
  ortho.big.3 <- merge(x = ortho.big.2, y = rbh, by = one.species, all = T)
  
  ## Checking interesting info
  
  # # How many genes are the same between Orthofinder orthologues and Blast RBH?
  # # for rows with single thaliana orthologues
  # sum(ortho.big.3$Arabidopsis_thaliana == ortho.big.3$RBH_B_id, na.rm = T) # 20448
  # # also for rows with multiple thaliana orthologues
  # # it might not work for rows with multiple target species RBH (the identical ones)
  # sum(mapply(FUN = grepl, pattern = ortho.big.3$RBH_B_id, x = ortho.big.3$Arabidopsis_thaliana), na.rm = T) # 20879
  # 
  # # Where Orthofinder orthologues and blast RBH differ?
  # sum(!mapply(FUN = grepl, pattern = ortho.big.3$RBH_B_id, x = ortho.big.3$Arabidopsis_thaliana), na.rm = T) # 426
  # View(ortho.big.3[!mapply(FUN = grepl, pattern = ortho.big.3$RBH_B_id, x = ortho.big.3$Arabidopsis_thaliana) &
  #               !is.na(ortho.big.3$Arabidopsis_thaliana) & 
  #               !is.na(ortho.big.3$RBH_B_id), ])
  
  
  ### Merging with simple Blast results
  
  t.print("Merging with simple Blast results for", one.species)
  
  # read Blast results
  
  # slow version
  # blast1 <- read.table(file = gzfile(paste0("blast_results/", one.species, "_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv.gz")), sep = "\t", header = F, comment.char = "#")
 
  # fast version
  # First pre-process the file in bash (remove comment lines using zgrep) and then read using fread.
  # fread doesn't support "comment.char".
  # I use linux in Windows (WSL), so I need the "wsl" in the cmd.
  blast1 <- fread(# file = (paste0("blast_results/", one.species, "_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv.gz")), 
                                cmd = paste0("wsl zgrep -v '^#' blast_results/", one.species, "_vs_Arabidopsis_thaliana_all_Vs_all_eval_0.1.tsv.gz"),
                                sep = "\t", header = F, 
                                # comment.char = "#"
                                )
  
  colnames(blast1) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qlen", "slen")
  # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, % query coverage per hsp, query length, subject length
  
  ## choose the best hits for each query (target species) gene
  
  t.print("choose the best hits for each query (target species) gene")
  
  ## fast version using data.table
  
  # get just the hits with highest bitscore for each query
  # if there are multiple such hits, collapse them to one line
  fn.max.bitscore <- function(x) {
    x[x$bitscore == max(x$bitscore), ][, lapply(X = .SD, FUN = function(x) {paste(unique(x), collapse = ", ")})]
  }
  
  blast.single.2 <- blast1[, fn.max.bitscore(.SD), by = qaccver]
  
  colnames(blast.single.2)[1] <- one.species
  colnames(blast.single.2)[-1] <- paste0("BLAST_", colnames(blast.single.2)[-1])
  
  # merging to large table
  ortho.big.4 <- merge(x = ortho.big.3, y = blast.single.2, by = one.species, all = T)
  
  
  
  ### Final selection of orthologues
  
  t.print("Final selection of orthologues for", one.species)
  
  ortho.big.5 <- ortho.big.4
  
  ## column with final orthologues
  
  # 1. use orthologues from orthologues table
  ortho.big.5$Arabidopsis_thaliana <- ortho.big.5$OL_Arabidopsis_thaliana
  # call them "orthologue"
  ortho.big.5$homologue_type <- NA
  ortho.big.5$homologue_type[!is.na(ortho.big.5$OL_Arabidopsis_thaliana)] <- "orthologue"
  
  # How many genes have orthologues?
  stats.ortho[1, "legend"] <- "A: OrthoFinder orthologues"
  stats.ortho[1, one.species] <- sum(!is.na(ortho.big.5$OL_Arabidopsis_thaliana))
  
  # 2. if missing, use Blast RBH
  ortho.big.5$Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana)] <- ortho.big.5$RBH_B_id[is.na(ortho.big.5$Arabidopsis_thaliana)] 
  
  # How many genes have orthologues?
  stats.ortho[2, "legend"] <- "B: BLAST Reciprocal Best Hits"
  stats.ortho[2, one.species] <- sum(!is.na(ortho.big.5$RBH_B_id))
  
  stats.ortho[3, "legend"] <- "A + B"
  stats.ortho[3, one.species] <- sum(!is.na(ortho.big.5$Arabidopsis_thaliana))
  
  # 3. if missing, use genes from N0 hierarchical orthogroups
  ortho.big.5$Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana)] <- ortho.big.5$N0_Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana)] 
  
  # How many genes have orthologues?
  stats.ortho[4, "legend"] <- "C: Shared N0 hierarchical orthogroups"
  stats.ortho[4, one.species] <- sum(!is.na(ortho.big.5$N0_Arabidopsis_thaliana))
  
  stats.ortho[5, "legend"] <- "A + B + C"
  stats.ortho[5, one.species] <- sum(!is.na(ortho.big.5$Arabidopsis_thaliana))
  
  # 4. if missing, use genes from broad orthogroups
  ortho.big.5$Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana)] <- ortho.big.5$OG_Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana)] 
  
  # How many genes have orthologues?
  stats.ortho[6, "legend"] <- "D: Shared Orthogroups"
  stats.ortho[6, one.species] <- sum(!is.na(ortho.big.5$OG_Arabidopsis_thaliana))
  
  stats.ortho[7, "legend"] <- "A + B + C + D"
  stats.ortho[7, one.species] <- sum(!is.na(ortho.big.5$Arabidopsis_thaliana))
  
  # 5. if missing, use genes hits from simple Blast (but only those with BLAST_pident > 40 & BLAST_qcovhsp > 50)
  ortho.big.5$Arabidopsis_thaliana[is.na(ortho.big.5$Arabidopsis_thaliana) & 
                                                                    ortho.big.5$BLAST_pident > 40 & 
                                                                    ortho.big.5$BLAST_qcovhsp > 50 &
                                                                    !is.na(ortho.big.5$BLAST_saccver)] <- 
    ortho.big.5$BLAST_saccver[is.na(ortho.big.5$Arabidopsis_thaliana) &
                                                            ortho.big.5$BLAST_pident > 40 & 
                                                            ortho.big.5$BLAST_qcovhsp > 50 &
                                                            !is.na(ortho.big.5$BLAST_saccver)] 
  
  # call other homologues than Orthofinder orthologue table orthologues as "other"
  ortho.big.5$homologue_type[is.na(ortho.big.5$homologue_type) & !is.na(ortho.big.5$Arabidopsis_thaliana)] <- "other"
  
  # How many genes have orthologues?
  stats.ortho[8, "legend"] <- "E: BLAST hits (pident > 40, qcovhsp > 50)"
  stats.ortho[8, one.species] <- sum(!is.na(ortho.big.5$BLAST_saccver) & 
                                       ortho.big.5$BLAST_pident > 40 & 
                                       ortho.big.5$BLAST_qcovhsp > 50)
  
  stats.ortho[9, "legend"] <- "A + B + C + D + E"
  stats.ortho[9, one.species] <- sum(!is.na(ortho.big.5$Arabidopsis_thaliana))
  
  
  ### Choosing the best orthologues for column with single orthologues
  
  ortho.big.5$single_Arabidopsis_thaliana <- NA
  
  # transfering the orthologues that are already single
  ortho.big.5$single_Arabidopsis_thaliana[!grepl(pattern = ", ", x = ortho.big.5$Arabidopsis_thaliana)] <- ortho.big.5$Arabidopsis_thaliana[!grepl(pattern = ", ", x = ortho.big.5$Arabidopsis_thaliana)]
  
  # Loop for choosing best orthologue/homologue if there are multiple
  
  t.print("loop for choosing best orthologue/homologue if there are multiple")
  
  for(i in grep(pattern = ", ", x = ortho.big.5$Arabidopsis_thaliana)) {
    gene.row <- ortho.big.5[i, ]
    split.ortho <- unlist(strsplit(x = gene.row$Arabidopsis_thaliana, split = ", "))
    # If there is one and only RBH and it is among orthologues, take the RBH
    if(!is.na(gene.row$RBH_B_id) & !grepl(pattern = ", ", x = gene.row$RBH_B_id) & gene.row$RBH_B_id %in% split.ortho) {
        # if(!gene.row$RBH_B_id %in% split.ortho) print(paste(gene.row$Arabidopsis_arenosa, ": RBH is not among orthologues"))
        s.ortho <- gene.row$RBH_B_id
    } else {
        # Else find the results from simple BLAST for this gene and find the orthologues among the hits
        sub.blast <- blast1[blast1$qaccver == as.character(gene.row[1, c(one.species), with = F]) & 
                                                    blast1$saccver %in% split.ortho, ]
        # If the orthologues are missing among BLAST hits, print warning.
        if(nrow(sub.blast) == 0) {
            print(paste0(gene.row[, one.species, with = F], ": Missing orthologue among BLAST hits"))
            s.orth <- NA
        } else {
            # Else take the hit with the highest bitscore. If several have the highest, it will take the first.
            s.ortho <- sub.blast$saccver[sub.blast$bitscore == max(sub.blast$bitscore)][1]
        }
    }
    ortho.big.5$single_Arabidopsis_thaliana[i] <- s.ortho
  }
  
  # how many genes without assigned orthologues do we have?
  stats.ortho[10, "legend"] <- "F: Genes without orthologues/homologues"
  stats.ortho[10, one.species] <- sum(is.na(ortho.big.5$Arabidopsis_thaliana))
  
  # Total number of genes
  stats.ortho[11, "legend"] <- "Genes total"
  stats.ortho[11, one.species] <- nrow(ortho.big.5)
  
  # Genes with "other" homologues
  stats.ortho[12, "legend"] <- "Genes with other homologues (B + C + D + E - A)"
  stats.ortho[12, one.species] <- sum(ortho.big.5$homologue_type == "other", na.rm = T)
  
  # how many genes multiple orthologues do we have?
  stats.ortho[13, "legend"] <- "Genes with multiple orthologues/homologues"
  stats.ortho[13, one.species] <- sum(grepl(pattern = ", ", x = ortho.big.5$Arabidopsis_thaliana))
  
  # reordering columns
  colnames(ortho.big.5)
  ortho.big.6 <- ortho.big.5[, c(1:2, 30, 32, 31, 3:29)]
  colnames(ortho.big.6)
  

  
  ### Exporting results
  
  t.print("Exporting tables for", one.species)
  
  dir.create(path = "supplemented_orthologues", recursive = T, showWarnings = F)
  write.table(x = ortho.big.6, file = paste0("supplemented_orthologues/", one.species, "__v__Arabidopsis_thaliana_supplemented.tsv"), sep = "\t", row.names = F)
  
  t.print(one.species, "done.")
  
  # Stop saving output to log file
  sink(type = "message")
  sink()
  close(log.file)
}
### End of the big loop


### Exporting table with statistics
dir.create(path = "figs_and_stats", recursive = T, showWarnings = F)
write.table(x = stats.ortho, file = "figs_and_stats/supplementing_orthologues_stats.tsv", sep = "\t", row.names = F)
```

# Explanation of columns in the output table

Name of output tables (in
[supplemented_orthologues](supplemented_orthologues/) directory):
`Species_name__v__Arabidopsis_thaliana_supplemented.tsv`

Columns:

- **Species_name**: Gene IDs of the investigated species.
- **prot_length**: Length of protein sequence of the longest/primary
  isoform used for OrthoFinder.
- **Arabidopsis_thaliana**: Selected *A. thaliana* orthologues according
  to these rules:  
  1. Use orthologues from orthologues table.  
  2. If missing, use Blast RBH.  
  3. If missing, use genes from N0 hierarchical orthogroups.  
  4. If missing, use genes from broad orthogroups.  
  5. If missing, use genes hits from simple Blast (but only those with
  BLAST_pident \> 40 & BLAST_qcovhsp \> 50)  
  BLAST_pident: Percentage of identical matches (in local alignment)  
  BLAST_qcovhsp: Query coverage per HSP (%)  
  HSP = High-scoring Segment Pair (local alignment with no gaps)  
- **single_Arabidopsis_thaliana**: Selected *A. thaliana* orthologue
  (only one) according to these rules:  
  1. If there is only one *A. thaliana* gene in column
  `Arabidopsis_thaliana`, use that one.  
  2. Else, if there is one and only RBH and it is among genes in column
  `Arabidopsis_thaliana`, take the RBH.  
  3. Else, take the gene from the column `Arabidopsis_thaliana` that had
  the highest bitscore in the simple BLAST results. If several have the
  highest bitscore, take the first.  
  4. If the genes in column `Arabidopsis_thaliana` are not among BLAST
  hits, take nothing.  
  **Recommendation**: It is appealing to use the “single” orthologues
  (`single_Arabidopsis_thaliana`), because it makes things easier.
  However, I recommend to do that only in cases when it is absolutely
  necessary and otherwise use the column where there are multiple
  orthologues in some cases (`Arabidopsis_thaliana`). The real orthology
  is not always one-to-one due to multiplications of genes in some
  species (there is a good explanation in [OrthoFinder
  GitHub](https://github.com/davidemms/OrthoFinder)). The column
  `Arabidopsis_thaliana` should better reflect the real biology.
- **homologue_type**: Type of selected *A. thaliana* orthologues in
  `Arabidopsis_thaliana` column.  
  `orthologue` - orthologue from OrthoFinder orthologues table (method
  1)  
  `other` - homologue assigned by other method (2-5)  
- **OL_Arabidopsis_thaliana**: *A. thaliana* orthologues from
  OrthoFinder orthologues table.
- **N0_HOG**: N0 orthogroup ID from OrthoFinder N0 hierarchical
  orthogroups table.
- **N0_Arabidopsis_thaliana**: *A. thaliana* orthologues/homologues from
  OrthoFinder N0 hierarchical orthogroups table.
- **OG_Orthogroup**: Orthogroup ID from OrthoFinder orthogroups table.
- **OG_Arabidopsis_thaliana**: *A. thaliana* orthologues/homologues from
  OrthoFinder orthogroups table.
- **RBH\_\[…\]**: Results of Reciprocal best hit (RBH) BLAST.  
  `RBH_B_id` - ID from species B (*A. thaliana*)  
  `RBH_A_length` - Length of protein sequence of the investigated (A)
  species (unfortunately including `*` in the end)  
  `RBH_B_length` - Length of *A. thaliana* (species B) hit protein
  sequence (unfortunately including `*` in the end)  
  `RBH_A_qcovhsp` - Percentage of sequence A covered by alignment
  (investigated species)  
  `RBH_B_qcovhsp` - Percentage of sequence B covered by alignment (*A.
  thaliana*)  
  `RBH_length` - HSP alignment length  
  `RBH_pident` - HSP percentage identity  
  `RBH_bitscore` - HSP bitscore  
- **BLAST\_\[…\]**: Results for highest-scoring hit of simple BLAST
  (investigated species against *A. thaliana* protein data  
  `BLAST_qaccver` - Query (investigated species) accession and version  
  `BLAST_saccver` - Subject (*A. thaliana*) accession and version  
  `BLAST_pident` - Percentage of identical matches  
  `BLAST_length` - Alignment length  
  `BLAST_mismatch` - Number of mismatches  
  `BLAST_gapopen` - Number of gap openings  
  `BLAST_qstart` - Start of alignment in query  
  `BLAST_qend` - End of alignment in query  
  `BLAST_sstart` - Start of alignment in subject (database hit)  
  `BLAST_send` - End of alignment in subject (database hit)  
  `BLAST_evalue` - Expectation value (E-value)  
  `BLAST_bitscore` - Bit score  
  `BLAST_qcovhsp` - Query coverage per HSP (alignment) \[%\]  
  `BLAST_qlen` - Query sequence length  
  `BLAST_slen` - Subject sequence length

# Improving table with stats

``` r
setwd("D:/!ecolgen/Brassicaceae_orthology/brassicaceae_3/")
if(!exists("old.par")) old.par<-par(no.readonly = T)

stats.ortho <- read.table(file = "figs_and_stats/supplementing_orthologues_stats.tsv", sep = "\t", head = T)

# stats.ortho.2 <- t(stats.ortho)

# put total on top
stats.ortho_2 <- stats.ortho[c(11, 1:10, 12, 13), ]
stats.ortho_3 <- stats.ortho_2
# calculate percent
for(i in 2:ncol(stats.ortho_2)) {
    stats.ortho_3[, i] <- paste0(formatC(x = stats.ortho_2[, i], big.mark = " "), " (", round(stats.ortho_2[, i]/stats.ortho_2[1, i], 2)*100, "%)")
}

write.table(x = stats.ortho_3, file = "figs_and_stats/supplementing_orthologues_stats_2.tsv", sep = "\t", row.names = F)
```

# End
