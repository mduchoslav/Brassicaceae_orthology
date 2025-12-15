SW installation for Annotation of *Alyssum gmelinii* genome
================
Miloš Duchoslav
2025-02

- [OrthoFinder](#orthofinder)
- [R](#r)
- [BLAST](#blast)
- [Reciprocal best hit (RBH) BLAST](#reciprocal-best-hit-rbh-blast)
- [Interproscan](#interproscan)

## OrthoFinder

version 2.5.4

<https://github.com/davidemms/OrthoFinder>

## R

Version: 4.3.3

``` r
R.Version()
```

## BLAST

I used BLAST installed as module on Metacentrum.

``` sh
# load module
module load blast-plus/2.16.0-gcc-10.2.1-bgzrrrz

# version
blastp -version
# blastp: 2.16.0+
#  Package: blast 2.16.0, build Nov  6 2024 15:13:03
```

## Reciprocal best hit (RBH) BLAST

Script from
<https://github.com/peterjc/galaxy_blast/blob/master/tools/blast_rbh>

``` bash
cd /storage/brno12-cerit/home/duchmil/SW

mkdir blast_rbh
cd blast_rbh/
# download script and description
wget https://raw.githubusercontent.com/peterjc/galaxy_blast/refs/heads/master/tools/blast_rbh/blast_rbh.py
wget https://raw.githubusercontent.com/peterjc/galaxy_blast/refs/heads/master/tools/blast_rbh/best_hits.py
wget https://raw.githubusercontent.com/peterjc/galaxy_blast/refs/heads/master/tools/blast_rbh/blast_rbh.xml

# help
python3 blast_rbh.py -h
```

## Interproscan

<https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html>

``` sh
cd /storage/brno12-cerit/home/duchmil/SW/InterProScan/

wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.76-107.0-64-bit.tar.gz.md5
# Must return *interproscan-5.76-107.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

tar -pxvzf interproscan-5.76-107.0-64-bit.tar.gz
# where:
#     p = preserve the file permissions
#     x = extract files from an archive
#     v = verbosely list the files processed
#     z = filter the archive through gzip
#     f = use archive file

cd interproscan-5.76-107.0/

#Index hmm models
python3 setup.py -f interproscan.properties

# check
module load openjdk
module list
# openjdk/17.0.0_35-intel-19.0.4-udmzouu
./interproscan.sh

# tests
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
./interproscan.sh -i test_all_appl.fasta -f tsv
```
