SW installation for Annotation of *Alyssum gmelinii* genome
================
Miloš Duchoslav
2025-02

- [OrthoFinder](#orthofinder)
- [R](#r)
- [BLAST](#blast)
- [Reciprocal best hit (RBH) BLAST](#reciprocal-best-hit-rbh-blast)

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
