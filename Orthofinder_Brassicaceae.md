Orthofinder - Brassicaceae
================
Miloš Duchoslav
2025-02

- [Introduction](#introduction)
- [Orthofinder Brassicaceae 2](#orthofinder-brassicaceae-2)
  - [Species included in this Orthofinder
    run](#species-included-in-this-orthofinder-run)
- [Preparation of protein sequences - Brassicaceae
  1](#preparation-of-protein-sequences---brassicaceae-1)
  - [Species with the downloadable primary
    transcripts](#species-with-the-downloadable-primary-transcripts)
  - [Species anotated at NCBI](#species-anotated-at-ncbi)
  - [Species other than annotated by NCBI with available protein
    sequences](#species-other-than-annotated-by-ncbi-with-available-protein-sequences)
  - [Species without downloadable protein
    sequences](#species-without-downloadable-protein-sequences)
  - [Extraction of primary transcripts using Orthofinder
    script](#extraction-of-primary-transcripts-using-orthofinder-script)
- [Preparation of protein sequences - Brassicaceae
  2](#preparation-of-protein-sequences---brassicaceae-2)
  - [Copying files from Brassicaceae
    1](#copying-files-from-brassicaceae-1)
  - [Protein sequences from new
    genomes](#protein-sequences-from-new-genomes)
  - [Extraction of primary transcripts using Orthofinder
    script](#extraction-of-primary-transcripts-using-orthofinder-script-1)
- [Orthofinder run](#orthofinder-run)
  - [brassicaceae_2 run](#brassicaceae_2-run)
- [Results from Orthofinder](#results-from-orthofinder)
  - [Visualization of statistics from
    Orthofinder](#visualization-of-statistics-from-orthofinder)

# Introduction

This RMarkdown file (or its markdown version for GitHub) documents
running of OrthoFinder tool to get orthology relationships between some
Brassicaceae species.

This file includes:

1.  BASH code that I ran at MetaCentrum (Czech national grid
    infrastructure) running PBS scheduling system for batch jobs.
2.  R code that I ran locally.

### SW installation and versions

The SW installation instructions and versions of SW used is described in
[Installation_of_SW.md](Installation_of_SW.md).

# Orthofinder Brassicaceae 2

This is run *Brassicaceae 2*. Compared to run *Brassicaceae 1*, I added
two species, whose genomes we recently assembled:  
- *Cardamine glauca*  
- *Noccaea praecox*

## Species included in this Orthofinder run

1.  *Alyssum gmelinii*
    - Celestini, Sonia, Miloš Duchoslav, Mahnaz Nezamivand-Chegini, Jörn
      Gerchen, Gabriela Šrámková, Raúl Wijfjes, Anna Krejčová, et
      al. “Genomic Basis of Adaptation to Serpentine Soil in Two Alyssum
      Species Shows Convergence with Arabidopsis across 20 Million Years
      of Divergence.” bioRxiv, February 28, 2025.
      <https://doi.org/10.1101/2025.02.27.640498>.
2.  *Arabidopsis arenosa*
    - Bramsiepe, Jonathan, Anders K. Krabberød, Katrine N. Bjerkan,
      Renate M. Alling, Ida M. Johannessen, Karina S. Hornslien,
      Jason R. Miller, Anne K. Brysting, and Paul E. Grini. “Structural
      Evidence for MADS-Box Type I Family Expansion Seen in New
      Assemblies of Arabidopsis Arenosa and A. Lyrata.” The Plant
      Journal 116, no. 3 (2023): 942–61.
      <https://doi.org/10.1111/tpj.16401>.
3.  *Arabidopsis lyrata* NCBI
    - Genome from Hu, Tina T., Pedro Pattyn, Erica G. Bakker, Jun Cao,
      Jan-Fang Cheng, Richard M. Clark, Noah Fahlgren, et al. “The
      Arabidopsis Lyrata Genome Sequence and the Basis of Rapid Genome
      Size Change.” Nature Genetics 43, no. 5 (May 2011): 476–81.
      <https://doi.org/10.1038/ng.807>.
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000004255.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004255.2/)
4.  *Arabidopsis lyrata* Rawat
    - Genome from Hu, Tina T., Pedro Pattyn, Erica G. Bakker, Jun Cao,
      Jan-Fang Cheng, Richard M. Clark, Noah Fahlgren, et al. “The
      Arabidopsis Lyrata Genome Sequence and the Basis of Rapid Genome
      Size Change.” Nature Genetics 43, no. 5 (May 2011): 476–81.
      <https://doi.org/10.1038/ng.807>.
    - Annotation from Rawat, Vimal, Ahmed Abdelsamad, Björn Pietzenuk,
      Danelle K. Seymour, Daniel Koenig, Detlef Weigel, Ales Pecinka,
      and Korbinian Schneeberger. “Improving the Annotation of
      Arabidopsis Lyrata Using RNA-Seq Data.” PLOS ONE 10, no. 9
      (September 18, 2015): e0137391.
      <https://doi.org/10.1371/journal.pone.0137391>.
5.  *Arabidopsis thaliana*
    - Araport11 protein sequences (version 2022-09-14) downloaded from
      [arabidopsis.org](https://www.arabidopsis.org/download/file?path=Proteins/Araport11_protein_lists/archived/Araport11_pep_20220914_representative_gene_model.gz)
6.  *Arabis alpina*
    - Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig,
      Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria
      Willing, et al. “Improving and Correcting the Contiguity of
      Long-Read Genome Assemblies of Three Plant Species Using Optical
      Mapping and Chromosome Conformation Capture Data.” Genome Research
      27, no. 5 (May 1, 2017): 778–86.
      <https://doi.org/10.1101/gr.213652.116>.
    - Data: <http://www.arabis-alpina.org/refseq.html>, I used the
      version 5.1 of the genome (later than Jiao et al. 2017)
7.  *Brassica oleracea*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000695525.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000695525.1/)
8.  *Brassica rapa*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000309985.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000309985.2/)
9.  *Camelina sativa*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000633955.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633955.1/)
10. *Capsella rubella*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000375325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000375325.1/)
11. *Cardamine glauca*
    - Nezamivand-Chegini Mahnaz, Duchoslav Miloš, …, Kolar Filip 2025
      (to be published)
12. *Cardamine hirsuta*
    - Gan, Xiangchao, Angela Hay, Michiel Kwantes, Georg Haberer, Asis
      Hallab, Raffaele Dello Ioio, Hugo Hofhuis, et al. “The Cardamine
      Hirsuta Genome Offers Insight into the Evolution of Morphological
      Diversity.” Nature Plants 2, no. 11 (October 31, 2016): 1–7.
      <https://doi.org/10.1038/nplants.2016.167>.
    - Data: <http://chi.mpipz.mpg.de/assembly.html>
13. *Cochlearia excelsa*
    - Bray, Sian M., Tuomas Hämälä, Min Zhou, Silvia Busoms, Sina
      Fischer, Stuart D. Desjardins, Terezie Mandáková, et
      al. “Kinetochore and Ionomic Adaptation to Whole-Genome
      Duplication in Cochlearia Shows Evolutionary Convergence in Three
      Autopolyploids.” Cell Reports 43, no. 8 (August 27, 2024).
      <https://doi.org/10.1016/j.celrep.2024.114576>.
    - Data: <https://doi.org/10.5061/dryad.ncjsxkt1s>
14. *Conringia planisiliqua*
    - Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig,
      Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria
      Willing, et al. “Improving and Correcting the Contiguity of
      Long-Read Genome Assemblies of Three Plant Species Using Optical
      Mapping and Chromosome Conformation Capture Data.” Genome Research
      27, no. 5 (May 1, 2017): 778–86.
      <https://doi.org/10.1101/gr.213652.116>.
15. *Euclidium syriacum*
    - Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig,
      Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria
      Willing, et al. “Improving and Correcting the Contiguity of
      Long-Read Genome Assemblies of Three Plant Species Using Optical
      Mapping and Chromosome Conformation Capture Data.” Genome Research
      27, no. 5 (May 1, 2017): 778–86.
      <https://doi.org/10.1101/gr.213652.116>.
16. *Eutrema salsugineum*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000478725.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000478725.1/)
17. *Noccaea praecox*
    - Nezamivand-Chegini Mahnaz, Duchoslav Miloš, …, Kolar Filip 2025
      (to be published)
18. *Raphanus sativus*
    - Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
    - [GCF_000801105.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000801105.1/)
19. *Rorippa islandica*
    - <https://phytozome-next.jgi.doe.gov/info/Rislandica_v1_1>

# Preparation of protein sequences - Brassicaceae 1

Here I describe the preparation of the data for the first run, then in
the next section, I will add the data for the second run.

### Folders etc.

``` sh
cd /storage/brno12-cerit/home/duchmil/orthofinder
mkdir brassicaceae_1
cd brassicaceae_1/
mkdir original_protein_fasta
mkdir renamed_protein_fasta_ncbi
mkdir renamed_protein_fasta_other
mkdir primary_transcripts
mkdir genomes_for_extracting_proteins
```

## Species with the downloadable primary transcripts

``` sh
### Download of the protein sequences
cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_1/original_protein_fasta

### Species with the downloadable primary transcripts

## A. thaliana
# downloading A. thaliana "representative gene model" from TAIR (Araport11)
wget https://www.arabidopsis.org/download_files/Proteins/Araport11_protein_lists/Araport11_pep_20220914_representative_gene_model.gz
gunzip Araport11_pep_20220914_representative_gene_model.gz
# the file is encoded in UTF-8
file -bi Araport11_pep_20220914_representative_gene_model
# There are sometimes strange characters and the OrthoFinder has with it problem. Thus, I will convert it to us-ascii, with removal of invalid characters (option -c).
iconv -c -f UTF-8 -t US-ASCII Araport11_pep_20220914_representative_gene_model -o Araport11_pep_20220914_representative_gene_model_cor.fasta
# renaming
cp Araport11_pep_20220914_representative_gene_model_cor.fasta ../primary_transcripts/Arabidopsis_thaliana.fasta

## Rorippa islandica
# https://phytozome-next.jgi.doe.gov/info/Rislandica_v1_1
# This curl command from JGI Data Portal is probably only short-term working:
curl --cookie jgi_session=/api/sessions/51ebfe2aa63f9400f9369bd81038003c --output R_islandica.zip -d "{\"ids\":{\"Phytozome-473\":[\"5f4058459a211ae42a1a2c03\",\"5f4058459a211ae42a1a2c01\"]}}" -H "Content-Type: application/json" https://files.jgi.doe.gov/filedownload/
unzip R_islandica.zip -d R_islandica
gunzip R_islandica/Phytozome/PhytozomeV13/Rislandica/v1.1/annotation/*.gz
cp R_islandica/Phytozome/PhytozomeV13/Rislandica/v1.1/annotation/Rislandica_473_v1.1.protein_primaryTranscriptOnly.fa ../primary_transcripts/Rorippa_islandica.fasta
```

## Species anotated at NCBI

``` sh
### Species anotated at NCBI
cd ~/orthofinder/brassicaceae_1/original_protein_fasta/

# loop for downloading of protein fastas and GFFs
while IFS=$'\t' read -r species_long_name species accession
do
# download from NCBI API
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"$accession"/download?include_annotation_type=PROT_FASTA,GENOME_GFF&filename="$accession".zip" -H "Accept: application/zip"
unzip $accession".zip" -d $accession
cp $accession"/ncbi_dataset/data/"$accession"/protein.faa" "../renamed_protein_fasta_ncbi/"$species".fasta"
cp $accession"/ncbi_dataset/data/"$accession"/genomic.gff" "../renamed_protein_fasta_ncbi/"$species".gff"
done < ../NCBI_annotated_genomes.txt # file with accession and species names to use


cd ../renamed_protein_fasta_ncbi/
# Rscript --verbose ../NCBI_longest_protein_executable.r B_oleracea.fasta B_oleracea.gff

# load R
module load r/4.1.3-gcc-10.2.1-6xt26dl

# loop for extracting the longest transcripts using my R script
species_names=$(cut -f 2 ../NCBI_annotated_genomes.txt)
echo $species_names
for species in $species_names
do
Rscript --verbose ../NCBI_longest_protein_executable.r $species".fasta" $species".gff"
done

# copy the files to the folder with all primary transcripts
cp longest_transcripts/*.fasta ../primary_transcripts/
```

## Species other than annotated by NCBI with available protein sequences

``` sh
### Species with available protein sequences other than annotated by NCBI
# Orthofinder script for extracting primary transcripts will be used for these.

cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_1/original_protein_fasta/

## A. arenosa
# https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_026151155.1/
# Protein sequences previously extracted from genome sequence and annotation by AGAT
cp ~/arenosa_genome_grini/arenosa_grini_proteins.fasta .
# renaming
cp arenosa_grini_proteins.fasta ../renamed_protein_fasta_other/Arabidopsis_arenosa.fasta

## A. lyrata - Rawat annotation
# Rawat et al., “Improving the Annotation of Arabidopsis Lyrata Using RNA-Seq Data.”
# <https://doi.org/10.1371/journal.pone.0137391>
# Release date: 2015
# Protein sequences previously extracted from genome sequence and modified annotation by AGAT.
# Copied from local place: d:\!ecolgen\resources\translation of gff to proteins\output\LyV2_proteins_mod_parents.fasta
# renaming
cp LyV2_proteins_mod_parents.fasta ../renamed_protein_fasta_other/Arabidopsis_lyrata_Rawat.fasta

## C. hirsuta
# Gan et al. “The Cardamine Hirsuta Genome Offers Insight into the Evolution of Morphological Diversity.” Nature Plants 2, no. 11 (October 31, 2016): 1–7. https://doi.org/10.1038/nplants.2016.167.
# http://chi.mpipz.mpg.de/assembly.html
wget http://chi.mpipz.mpg.de/download/annotations/carhr38.aa.fa

# Adding gene ID to fasta headers (It is needed for Orthofinder script to get longest isoforms.)
# (awk by ChatGPT)
awk '/^>/ {print $0 " gene=" gensub(/^>([A-Za-z0-9]+)\.[0-9]+/, "\\1", "g"); next} {print}' carhr38.aa.fa | head -n 100
awk '/^>/ {print $0 " gene=" gensub(/^>([A-Za-z0-9]+)\.[0-9]+/, "\\1", "g"); next} {print}' carhr38.aa.fa > Cardamine_hirsuta_with_gene_IDs.fasta
cp Cardamine_hirsuta_with_gene_IDs.fasta ../renamed_protein_fasta_other/Cardamine_hirsuta.fasta
# wget http://chi.mpipz.mpg.de/download/annotations/carhr38.gff
# cp carhr38.gff ../renamed_protein_fasta_other/Cardamine_hirsuta.gff

## Cochlearia excelsa
# Bray, Sian M., Tuomas Hämälä, Min Zhou, Silvia Busoms, Sina Fischer, Stuart D. Desjardins, Terezie Mandáková, et al. “Kinetochore and Ionomic Adaptation to Whole-Genome Duplication in Cochlearia Shows Evolutionary Convergence in Three Autopolyploids.” Cell Reports 43, no. 8 (August 27, 2024). https://doi.org/10.1016/j.celrep.2024.114576.
# Data: https://doi.org/10.5061/dryad.ncjsxkt1s

# I was not able to download from Dryad using wget, I had to download the files using web browser.
gunzip C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta.gz
grep -c '>' C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta # 54424

# It seems like there are only primary transcripts, I probably don't need to extract them.
grep -P -c "\tgene\t" C_excelsa_V5_braker2_wRseq.gff3 # 54424
grep -P -c "\tmRNA\t" C_excelsa_V5_braker2_wRseq.gff3 # 56280

# Adding gene ID to fasta headers (It is needed for Orthofinder script to get longest isoforms.)
# (awk by ChatGPT)
awk '/^>/ {print $0 " gene=" gensub(/^>([A-Za-z0-9]+)\.t[0-9]+/, "\\1", "g"); next} {print}' C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta | head -n 100
awk '/^>/ {print $0 " gene=" gensub(/^>([A-Za-z0-9]+)\.t[0-9]+/, "\\1", "g"); next} {print}' C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta > Cochlearia_excelsa_with_gene_IDs.fasta
cp Cochlearia_excelsa_with_gene_IDs.fasta ../renamed_protein_fasta_other/Cochlearia_excelsa.fasta





## Alyssum gmelinii (A. montanum group)
# Genome assembly and annotation:
# Sonia Celestini, Miloš Duchoslav, Mahnaz Nezamivand-Chegini, Jörn Gerchen, Gabriela Šrámková, Raúl Wijfjes, Anna Krejčová, Nevena Kuzmanović, Stanislav Španiel, Korbinian Schneeberger, Levi Yant, Filip Kolář (2025): Genomic basis of adaptation to serpentine soil in two *Alyssum* species shows convergence with *Arabidopsis* across 20 million years of divergence. Preprint at bioRxiv: <https://doi.org/10.1101/2025.02.27.640498>
cp /storage/brno12-cerit/home/duchmil/annotations/alyssum_2024_Mahnaz_assembly/annot_processing/alyssum_v2_annotation_proteins.fasta .
cp alyssum_v2_annotation_proteins.fasta ../renamed_protein_fasta_other/Alyssum_gmelinii.fasta
```

## Species without downloadable protein sequences

### Downloading genomes and annotations

``` sh
### Species without downloadable protein sequences
cd ~/orthofinder/brassicaceae_1/genomes_for_extracting_proteins

## Arabis alpina
# Jiao et al. 2017 (https://doi.org/10.1101%2Fgr.213652.116)
# http://www.arabis-alpina.org/refseq.html, I used the version 5.1 of the genome (later than Jiao et al. 2017)
species=Arabis_alpina
mkdir $species
cd $species
wget http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz
wget http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/A.alpina.modified.5.1.gff3.gz
# AGAT can read gzipped GFF, but not gzipped genome fasta
gunzip Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz
cd ..

# Download annotation of Conringia planisiliqua and Euclidium syriacum
# from supplementary files of Jiao et al. 2017 (https://doi.org/10.1101%2Fgr.213652.116)
wget https://genome.cshlp.org/content/suppl/2017/04/05/gr.213652.116.DC1/Supplemental_Data.zip
mv Supplemental_Data.zip Jiao_2017_Supplemental_Data.zip
unzip Jiao_2017_Supplemental_Data.zip
rm Arabis_alpina.MPIPZ.version_5.chr.all.liftOverV4.v3.gff3.gz

## Conringia planisiliqua
# Jiao et al. 2017 (https://doi.org/10.1101%2Fgr.213652.116)
# Genome annotation (GFF) is from supplementary data of Jiao et al. 2017.
accession=GCA_900108845.1
species=Conringia_planisiliqua
mkdir $species
cd $species
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"$accession"/download?include_annotation_type=GENOME_FASTA&filename="$accession".zip" -H "Accept: application/zip"
unzip $accession".zip" -d $accession
cp $accession"/ncbi_dataset/data/"$accession"/"*".fna" ./$species"_genome.fasta"
# Geting annotation from downloaded supplements
mv ../$species*.gff.gz .
# Unzipping genome annotation
gunzip $species*.gff.gz
cd ..

## Euclidium syriacum
# Jiao et al. 2017 (https://doi.org/10.1101%2Fgr.213652.116)
# Genome annotation (GFF) is from supplementary data of Jiao et al. 2017.
accession=GCA_900116095.1
species=Euclidium_syriacum
mkdir $species
cd $species
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"$accession"/download?include_annotation_type=GENOME_FASTA&filename="$accession".zip" -H "Accept: application/zip"
unzip $accession".zip" -d $accession
cp $accession"/ncbi_dataset/data/"$accession"/"*".fna" ./$species"_genome.fasta"
# Geting annotation from downloaded supplements
mv ../$species*.gff.gz .
# Unzipping genome annotation
gunzip $species*.gff.gz
cd ..


# There is a problem that the IDs (scaffold names) are different in annotation and in genome sequence. This I will solve locally in R. I downloaded the files.
```

### Changing of sequence names in fasta files with genomes of E. syriacum and C. planisiliqua

``` r
setwd("D:/!ecolgen/resources/orthofinder/brassicaceae_1")
```

``` r
### C. planisiliqua

## Changing sequence names in fasta genome file to match those in GFF annotation

# Fasta to read
arg.fasta <- "genomes_for_extracting_proteins/Conringia_planisiliqua/Conringia_planisiliqua_genome.fasta"
# Fasta to write
fasta.out <- "genomes_for_extracting_proteins/Conringia_planisiliqua/Conringia_planisiliqua_genome_seqnames_as_gff.fasta"

# reading fasta
library(seqinr)
genome <- read.fasta(file = arg.fasta, as.string = T, forceDNAtolower = F, strip.desc = T)

# extracting annotation of seqeunces
seq.annot <- sapply(X = genome, FUN = attr, which = "Annot")
# extracting scaffold name form annotations
seq.names <- gsub(pattern = "(^.*contig: )|(, whole.*$)", replacement = "", x = seq.annot)

# writing fasta file
write.fasta(sequences = genome, names = paste(seq.names, seq.annot), file.out = fasta.out, nbchar = 80, as.string = T)


## Check whether the GFF file seems to match the genome assembly

# read the gff file
gff.2 <- read.table(file = "genomes_for_extracting_proteins/Conringia_planisiliqua/Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff",
                    header = F, sep = "\t", comment.char = "#", quote = "" 
                    # , nrows = 1500
                    )
summary(gff.2)
head(gff.2)
levels(as.factor(gff.2$V3))

## Checking the length of scaffolds from GFF and fasta

# length of scaffolds from fasta
seq.length <- getLength(genome)

# maximal coordinates on each scaffold in GFF
seq.length.gff <- tapply(X = gff.2$V5, INDEX = gff.2$V1, FUN = max)

# comparison
seq.l.df <- cbind.data.frame(seq.names, seq.length, seq.length.gff[match(seq.names, names(seq.length.gff))])
seq.l.df.ordered <- seq.l.df[order(row.names(seq.l.df)), ]
seq.l.df.ordered
# Is the lenght of scaffold always higher then the maximal coordinate?
all(seq.l.df$seq.length > seq.l.df$`seq.length.gff[match(seq.names, names(seq.length.gff))]`, na.rm = T)



### E. syriacum

## Changing sequence names in fasta genome file to match those in GFF annotation

# Fasta to read
arg.fasta <- "genomes_for_extracting_proteins/Euclidium_syriacum/Euclidium_syriacum_genome.fasta"
# Fasta to write
fasta.out <- "genomes_for_extracting_proteins/Euclidium_syriacum/Euclidium_syriacum_genome_seqnames_as_gff.fasta"

# reading fasta
library(seqinr)
genome <- read.fasta(file = arg.fasta, as.string = T, forceDNAtolower = F, strip.desc = T)

# extracting annotation of seqeunces
seq.annot <- sapply(X = genome, FUN = attr, which = "Annot")
# extracting scaffold name form annotations
seq.names <- gsub(pattern = "(^.*contig: )|(, whole.*$)", replacement = "", x = seq.annot)

# writing fasta file
write.fasta(sequences = genome, names = paste(seq.names, seq.annot), file.out = fasta.out, nbchar = 80, as.string = T)




## Check whether the GFF file seems to match the genome assembly

# read the gff file
gff.2 <- read.table(file = "genomes_for_extracting_proteins/Euclidium_syriacum/Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff",
                    header = F, sep = "\t", comment.char = "#", quote = "" 
                    # , nrows = 1500
                    )
summary(gff.2)
head(gff.2)
levels(as.factor(gff.2$V3))
table(gff.2$V3)

## Checking the length of scaffolds from GFF and fasta

# length of scaffolds from fasta
seq.length <- getLength(genome)

# maximal coordinates on each scaffold in GFF
seq.length.gff <- tapply(X = gff.2$V5, INDEX = gff.2$V1, FUN = max)

# comparison
seq.l.df <- cbind.data.frame(seq.names, seq.length, seq.length.gff[match(seq.names, names(seq.length.gff))])
seq.l.df.ordered <- seq.l.df[order(row.names(seq.l.df)), ]
seq.l.df.ordered
# Is the lenght of scaffold always higher then the maximal coordinate?
all(seq.l.df$seq.length > seq.l.df$`seq.length.gff[match(seq.names, names(seq.length.gff))]`, na.rm = T)
```

The gff files seem to match the assemblies.

### Extraction of protein sequences from genome and annotation (A. alpina, E. syriacum and C. planisiliqua) using AGAT

``` bash
# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_1/genomes_for_extracting_proteins


cd Arabis_alpina/

zgrep -P -c "\ttransposable_element_gene\t" A.alpina.modified.5.1.gff3.gz # 0
zgrep -P -c "\tgene\t" A.alpina.modified.5.1.gff3.gz # 34212

# run the AGAT container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

## Protein sequences
# AGAT can read gzipped GFF, but not gzipped genome fasta
agat_sp_extract_sequences.pl -g A.alpina.modified.5.1.gff3.gz -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta -p -o Arabis_alpina_proteins.fasta

# counting the number of sequences
grep ">" Arabis_alpina_proteins.fasta -c
# 34212

# converting from folded fasta to unfolded fasta for better counting and checking for internal stop codons
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Arabis_alpina_proteins.fasta | grep \*[[:alpha:]] | wc -l
# There are 463 proteins with internal stop codons.

# checking for alternative transcripts
grep ".t2 " Arabis_alpina_proteins.fasta -c # 0
# It seems that there is always only one transcript per gene.

# copy protein sequences
# There are no alternative splicings so I can put it directly to primary_transcripts
cp Arabis_alpina_proteins.fasta ../../primary_transcripts/Arabis_alpina.fasta

cd ..


cd Conringia_planisiliqua

# The annotation contains also transposable_element_gene features with mRNA_TE_gene, exon and CDS subfeatures. 
# They have to be removed, otherwise the proteins encoded by these genes will be also extracted.
grep -P -c "\ttransposable_element_gene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff # 21515
grep -P -c "\tgene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff # 34766

# list of protein coding genes
grep -P "\tgene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff | cut -f 9 | sed "s/ID=//" | sed "s/;.*$//" | head
grep -P "\tgene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff | cut -f 9 | sed "s/ID=//" | sed "s/;.*$//" > Conringia_planisiliqua_protein_coding_genes.txt
wc -l Conringia_planisiliqua_protein_coding_genes.txt # 34766

## filter to keep only protein coding genes and remove TEs
agat_sp_filter_feature_from_keep_list.pl --gff Conringia_planisiliqua.MPIPZ_v1.annotation.TE.gff --keep_list Conringia_planisiliqua_protein_coding_genes.txt -o Conringia_planisiliqua.MPIPZ_v1.annotation.without_TE.gff

# exit AGAT container
exit

grep -P -c "\ttransposable_element_gene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.without_TE.gff # 0
grep -P -c "\tgene\t" Conringia_planisiliqua.MPIPZ_v1.annotation.without_TE.gff # 34766

# run the AGAT container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

## Extract protein sequences
agat_sp_extract_sequences.pl -g Conringia_planisiliqua.MPIPZ_v1.annotation.without_TE.gff -f Conringia_planisiliqua_genome_seqnames_as_gff.fasta -p -o Conringia_planisiliqua_proteins.fasta

# exit AGAT container
exit

# counting the number of sequences
grep ">" Conringia_planisiliqua_proteins.fasta -c
# 34766

# converting from folded fasta to unfolded fasta for better counting and checking for internal stop codons
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Conringia_planisiliqua_proteins.fasta | grep \*[[:alpha:]] | wc -l
# There are 17 proteins with internal stop codons.

# checking for alternative transcripts
grep ".t2 " Conringia_planisiliqua_proteins.fasta -c # 0
# It seems that there is always only one transcript per gene.

# copy protein sequences
# There are no alterbative splicings so I can put it directly to primary_transcripts
cp Conringia_planisiliqua_proteins.fasta ../../primary_transcripts/Conringia_planisiliqua.fasta

cd ..


cd Euclidium_syriacum

# The annotation contains also transposable_element_gene features with mRNA_TE_gene, exon and CDS subfeatures. 
# They have to be removed, otherwise the proteins encoded by these genes will be also extracted.
grep -P -c "\ttransposable_element_gene\t" Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff # 33000
grep -P -c "\tgene\t" Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff # 33001

# list of protein coding genes
grep -P "\tgene\t" Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff | cut -f 9 | sed "s/ID=//" | sed "s/;.*$//" | head
grep -P "\tgene\t" Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff | cut -f 9 | sed "s/ID=//" | sed "s/;.*$//" > Euclidium_syriacum_protein_coding_genes.txt
wc -l Euclidium_syriacum_protein_coding_genes.txt # 33001

# run the AGAT container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

## filter to keep only protein coding genes and remove TEs
agat_sp_filter_feature_from_keep_list.pl --gff Euclidium_syriacum.MPIPZ_v1.annotation.TE.gff --keep_list Euclidium_syriacum_protein_coding_genes.txt -o Euclidium_syriacum.MPIPZ_v1.annotation.without_TE.gff

## Protein sequences
# running command and at the same time putting the messages to log file
agat_sp_extract_sequences.pl -g Euclidium_syriacum.MPIPZ_v1.annotation.without_TE.gff -f Euclidium_syriacum_genome_seqnames_as_gff.fasta -p -o Euclidium_syriacum_proteins.fasta

# exit AGAT container
exit

# counting the number of sequences
grep ">" Euclidium_syriacum_proteins.fasta -c
# 33001

# converting from folded fasta to unfolded fasta for better counting and checking for internal stop codons
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Euclidium_syriacum_proteins.fasta | grep \*[[:alpha:]] | wc -l
# There are 2 proteins with internal stop codons.

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Euclidium_syriacum_proteins.fasta | grep \*[[:alpha:]] -b1

# checking for alternative transcripts
grep ".t2 " Euclidium_syriacum_proteins.fasta -c # 0
# It seems that there is always only one transcript per gene.

# copy protein sequences
# There are no alterbative splicings so I can put it directly to primary_transcripts
cp Euclidium_syriacum_proteins.fasta ../../primary_transcripts/Euclidium_syriacum.fasta
```

## Extraction of primary transcripts using Orthofinder script

``` bash
### Extraction of primary transcripts using Orthofinder script

# loading Python and the modules needed
# module load python
module add python36-modules-gcc

# test that Orthofinder is running
~/orthofinder/OrthoFinder_source/orthofinder.py -h

## taking only the longest transcript using Orthofinder script
cd ~/orthofinder/brassicaceae_1/renamed_protein_fasta_other/
for f in *.fasta ; do python ~/orthofinder/OrthoFinder_source/tools/primary_transcript.py $f ; done
# Note: This does not work properly for NCBI data.

# This did not remove any transcripts from Cochlearia_excelsa, however, there were probably already only primary transcripts (see before).


# Copying files
cp -v primary_transcripts/*.fasta ../primary_transcripts/
```

# Preparation of protein sequences - Brassicaceae 2

### Folders etc.

``` sh
cd /storage/brno12-cerit/home/duchmil/orthofinder
mkdir brassicaceae_2
cd brassicaceae_2/
mkdir original_protein_fasta
mkdir renamed_protein_fasta_ncbi
mkdir renamed_protein_fasta_other
mkdir primary_transcripts
mkdir genomes_for_extracting_proteins
```

## Copying files from Brassicaceae 1

``` sh
cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/primary_transcripts

cp -v /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_1/primary_transcripts/*.fasta .
```

## Protein sequences from new genomes

``` sh
cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/original_protein_fasta

## Cardamine glauca
# Genome assembly and annotation:
# Nezamivand-Chegini Mahnaz, Duchoslav Milos, ..., Kolar Filip 2025 (to be published)
cp -v /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2024_12/annot_processing/Cardamine_glauca_v1_annotation_proteins.fasta .
cp -v Cardamine_glauca_v1_annotation_proteins.fasta ../renamed_protein_fasta_other/Cardamine_glauca.fasta

## Noccaea praecox
# Genome assembly and annotation:
# Nezamivand-Chegini Mahnaz, Duchoslav Milos, ..., Kolar Filip 2025 (to be published)
cp -v /storage/brno12-cerit/home/duchmil/annotations/Noccaea_praecox_2024_12/annot_processing/Noccaea_praecox_v1_annotation_proteins.fasta .
cp -v Noccaea_praecox_v1_annotation_proteins.fasta ../renamed_protein_fasta_other/Noccaea_praecox.fasta
```

## Extraction of primary transcripts using Orthofinder script

``` bash
### Extraction of primary transcripts using Orthofinder script

# loading Python and the modules needed
# module load python
module add python36-modules-gcc

# test that Orthofinder is running
~/orthofinder/OrthoFinder_source/orthofinder.py -h

## taking only the longest transcript using Orthofinder script
cd ~/orthofinder/brassicaceae_2/renamed_protein_fasta_other/
for f in *.fasta ; do python ~/orthofinder/OrthoFinder_source/tools/primary_transcript.py $f ; done
# Note: This does not work properly for NCBI data.

# This did not remove any transcripts from Cochlearia_excelsa, however, there were probably already only primary transcripts (see before).


# Copying files
cp -v primary_transcripts/*.fasta ../primary_transcripts/
```

# Orthofinder run

## brassicaceae_2 run

18 species, 19 input files (two different annotations for
Arabidopsis_lyrata)  
Input fasta files:

Alyssum_gmelinii.fasta Arabidopsis_arenosa.fasta
Arabidopsis_lyrata_NCBI.fasta Arabidopsis_lyrata_Rawat.fasta
Arabidopsis_thaliana.fasta Arabis_alpina.fasta Brassica_oleracea.fasta
Brassica_rapa.fasta Camelina_sativa.fasta Capsella_rubella.fasta
Cardamine_glauca.fasta Cardamine_hirsuta.fasta Cochlearia_excelsa.fasta
Conringia_planisiliqua.fasta Euclidium_syriacum.fasta
Eutrema_salsugineum.fasta Noccaea_praecox.fasta Raphanus_sativus.fasta
Rorippa_islandica.fasta

``` bash
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

# Resources used: 7 h, 54% CPU, 11 GB memory.
```

``` bash
### brassicaceae_2 run
cd /storage/brno12-cerit/home/duchmil/orthofinder/brassicaceae_2/metacentrum_scripts/

# putting the job in the queue
qsub orthofinder_brassicaceae_2.bash

qstat -u duchmil
```

# Results from Orthofinder

## Visualization of statistics from Orthofinder

### Orthogroups - species overlap

``` r
setwd("D:/!ecolgen/resources/orthofinder/brassicaceae_2/")
old.par<-par(no.readonly = T)

spec.overlap <- read.table(file = "orthofinder_results/Results_brassicaceae_2/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv")

## heatmap with values
pdf ("R_analysis/Orthogroups_SpeciesOverlaps_heatmap.pdf", width=14, height=7, onefile = T)
par(mar = c(2, 12, 12, 2) + 0.1)
# input data
gdata <- spec.overlap
# turn it from similarity to distances
gdata.dist <- max(gdata) - gdata
# clustering
clust <- hclust(d = as.dist(gdata.dist), method = "single") # this method seemed best
# plot(clust)
# reordering according to clustering
gdata.clust <- gdata[clust$order, clust$order]
# reversing order
# gdata.rev <- gdata.clust[nrow(gdata.clust):1, ]
gdata.rev <- gdata.clust[, ncol(gdata.clust):1]
# image (simple heat map - heatmap has some strange coordinate system, so it is not useful)
image(t(as.matrix(gdata.rev)), xaxt = "n", yaxt = "n",
      main = "")
title(main = "Number of shared orthogroups", line = 7)
xaxs <- ((1:ncol(gdata.rev))-1)/(ncol(gdata.rev)-1)
yaxs <- ((1:nrow(gdata.rev))-1)/(nrow(gdata.rev)-1)
axis(side = 3, at = xaxs, labels = colnames(gdata.rev), las = 2, cex.axis = 0.9)
axis(side = 2, at = yaxs, labels = rownames(gdata.rev), tick = T, las = 2, cex.axis = 1)
# values
text(x = rep(x = xaxs, each = length(yaxs)), y = rep(x = yaxs, times = length(xaxs)), 
     labels = format(as.matrix(gdata.rev), big.mark = " "), 
     cex = 0.8)
par(old.par)
dev.off()



# clust <- hclust(d = as.dist(gdata.dist), method = "average")
# plot(clust)
# clust <- hclust(d = as.dist(gdata.dist), method = "ward.D")
# plot(clust)
# clust <- hclust(d = as.dist(gdata.dist), method = "ward.D2")
# plot(clust)
# 
# clust <- hclust(d = as.dist(gdata.dist), method = "complete")
# clust <- hclust(d = as.dist(gdata.dist), method = "mcquitty")
# clust <- hclust(d = as.dist(gdata.dist), method = "median")
# clust <- hclust(d = as.dist(gdata.dist), method = "centroid")
# clust$order
# plot(clust)
# colnames(gdata)[clust$order]
```

### Orthologues statistics

``` r
folder <- "orthofinder_results/Results_brassicaceae_2/Comparative_Genomics_Statistics/"
file.pattern <- "OrthologuesStats"
files <- list.files(folder, pattern = file.pattern)

pdf ("R_analysis/Comparative_Genomics_Statistics_heatmaps.pdf", width=14, height=7, onefile = T)
par(mar = c(2, 12, 12, 2) + 0.1)
for(i in 5:1) {
  gdata <- read.table(file = paste0(folder, files[i]))
  # reverse rows of the matrix to be shown in the same order as the original one on picture
  gdata.rev <- gdata[nrow(gdata):1, ]
  # image (simple heat map - heatmap has some strange coordinate system, so it is not useful)
  image(t(as.matrix(gdata.rev)), xaxt = "n", yaxt = "n",
        main = "")
  title(main = sub(pattern = "\\.tsv", replacement = "", x = files[i]), line = 7)
  xaxs <- ((1:ncol(gdata.rev))-1)/(ncol(gdata.rev)-1)
  yaxs <- ((1:nrow(gdata.rev))-1)/(nrow(gdata.rev)-1)
  axis(side = 3, at = xaxs, labels = colnames(gdata.rev), las = 2, cex.axis = 0.9)
  axis(side = 2, at = yaxs, labels = rownames(gdata.rev), tick = T, las = 2, cex.axis = 1)
  # values
  text(x = rep(x = xaxs, each = length(yaxs)), y = rep(x = yaxs, times = length(xaxs)), 
       labels = format(as.matrix(gdata.rev), big.mark = " "), 
       cex = 0.8)
  mtext(text = "Note: Species in rows against species in column.", side = 1, cex = 0.8)
}
par(old.par)
dev.off()
```

### Statistics per species

``` r
per.spec <- read.table(file = "orthofinder_results/Results_brassicaceae_2/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv", sep = "\t", quote = "", header = T, row.names = 1, nrows = 10)


## heatmap with values
pdf ("R_analysis/Statistics_PerSpecies_heatmap.pdf", width=14, height=7, onefile = T)
par(mar = c(2, 22, 12, 2) + 0.1)
# input data
gdata <- per.spec
# relative value per row for heatmap
gdata.rel <- gdata/rowMeans(gdata)
# rowMeans(gdata.rel)
# reversing order
gdata.rev.rel <- gdata.rel[nrow(gdata.rel):1, ]
gdata.rev <- gdata[nrow(gdata):1, ]
# gdata.rev <- gdata[, ncol(gdata):1]
# image (simple heat map - heatmap has some strange coordinate system, so it is not useful)
image(t(as.matrix(gdata.rev.rel)), xaxt = "n", yaxt = "n",
      main = "")
title(main = "Orthofinder_brassicaceae_2", line = 7)
xaxs <- ((1:ncol(gdata.rev))-1)/(ncol(gdata.rev)-1)
yaxs <- ((1:nrow(gdata.rev))-1)/(nrow(gdata.rev)-1)
axis(side = 3, at = xaxs, labels = colnames(gdata.rev), las = 2, cex.axis = 0.9)
axis(side = 2, at = yaxs, labels = rownames(gdata.rev), tick = T, las = 2, cex.axis = 1)
# values
text(x = rep(x = xaxs, each = length(yaxs)), y = rep(x = yaxs, times = length(xaxs)), 
     labels = format(round(as.matrix(gdata.rev)), big.mark = " "), 
     cex = 0.8)
par(old.par)
dev.off()
```
