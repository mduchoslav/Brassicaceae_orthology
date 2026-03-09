[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18924985.svg)](https://doi.org/10.5281/zenodo.18924985)

Orthology and functional annotation of some Brassicaceae genomes
================
Miloš Duchoslav
2025

- [Introduction](#Introduction)
- [Version history](#Version-history)
- [Species included in this version](#Species-included-in-this-version)
- [Scripts and detailed description of methods](#Scripts-and-detailed-description-of-methods)
- [Intermediate results](#Intermediate-results)
- [Final results](#Final-results)

## Introduction

The purpose of this project is to obtain gene orthology relationships between several Brassicaceae species and good functional annotation of their genes. It will be used as a resource in other projects of our team ([Plant ecological genomics group of Filip Kolář](https://www.plantecologicalgenomics.cz/), Charles University, Prague, Czech Republic), but it could be useful also for others. I focus on species that we use in our projects, but I added also some other genomes with reasonably good assembly and annotation quality to improve orthology inference.

## Version history

[Version 3.0](https://github.com/mduchoslav/Brassicaceae_orthology/tree/3.0) (2025-12) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18924985.svg)](https://doi.org/10.5281/zenodo.18924985)
- Species added:
	- *Aethionema saxatile*
	- *Brassica napus*
	- *Erysimum linariifolium*
	- *Odontarrhena muralis*
- Species with changed assembly or annotation:
	- *Arabidopsis thaliana* (version of data is now Araport11 2025-04-11)
	- *Alyssum gmelinii* (assembly is now V3)
	- *Cardamine glauca* (annotation version is now 1.1)
	- *Noccaea praecox* (annotation version is now 1.1)
- Species removed
	- *Rorippa islandica* (because of restrictions on dataset usage)


[Version 2.0](https://github.com/mduchoslav/Brassicaceae_orthology/tree/2.0) (2025-03)
- Added 2 species, whose genomes we recently assembled:
	- *Cardamine glauca*  
	- *Noccaea praecox*

## Species included in this version

1. *Aethionema saxatile*
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (to be published)
	- Available at <https://github.com/mduchoslav/Genome_annotations_ultramafic_Brassicaceae>
2. *Alyssum gmelinii*
	- Assembly V3 (compared to V2 the assembly included HiC data)
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (unpublished)
3. *Arabidopsis arenosa*
	- Bramsiepe, Jonathan, Anders K. Krabberød, Katrine N. Bjerkan, Renate M. Alling, Ida M. Johannessen, Karina S. Hornslien, Jason R. Miller, Anne K. Brysting, and Paul E. Grini. “Structural Evidence for MADS-Box Type I Family Expansion Seen in New Assemblies of Arabidopsis Arenosa and A. Lyrata.” The Plant Journal 116, no. 3 (2023): 942–61. <https://doi.org/10.1111/tpj.16401>.
4. *Arabidopsis lyrata* NCBI
	- Genome from Hu, Tina T., Pedro Pattyn, Erica G. Bakker, Jun Cao, Jan-Fang Cheng, Richard M. Clark, Noah Fahlgren, et al. “The Arabidopsis Lyrata Genome Sequence and the Basis of Rapid Genome Size Change.” Nature Genetics 43, no. 5 (May 2011): 476–81. <https://doi.org/10.1038/ng.807>.
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000004255.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004255.2/)
5. *Arabidopsis lyrata* Rawat
	- Genome from Hu, Tina T., Pedro Pattyn, Erica G. Bakker, Jun Cao, Jan-Fang Cheng, Richard M. Clark, Noah Fahlgren, et al. “The Arabidopsis Lyrata Genome Sequence and the Basis of Rapid Genome Size Change.” Nature Genetics 43, no. 5 (May 2011): 476–81. <https://doi.org/10.1038/ng.807>.
	- Annotation from Rawat, Vimal, Ahmed Abdelsamad, Björn Pietzenuk, Danelle K. Seymour, Daniel Koenig, Detlef Weigel, Ales Pecinka, and Korbinian Schneeberger. “Improving the Annotation of Arabidopsis Lyrata Using RNA-Seq Data.” PLOS ONE 10, no. 9 (September 18, 2015): e0137391. <https://doi.org/10.1371/journal.pone.0137391>.
6. *Arabidopsis thaliana*
	- Araport11 protein sequences (version 2025-04-11) downloaded from [arabidopsis.org](https://www.arabidopsis.org/download/file?path=Proteins/Araport11_protein_lists/Araport11_pep_20250411_representative_gene_model.gz)
7. *Arabis alpina*
	- Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig, Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria Willing, et al. “Improving and Correcting the Contiguity of Long-Read Genome Assemblies of Three Plant Species Using Optical Mapping and Chromosome Conformation Capture Data.” Genome Research 27, no. 5 (May 1, 2017): 778–86. <https://doi.org/10.1101/gr.213652.116>.
	- Data: <http://www.arabis-alpina.org/refseq.html>, I used the version 5.1 of the genome (later than Jiao et al. 2017)
8. *Brasica napus*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_020379485.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020379485.1/)
9. *Brassica oleracea*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000695525.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000695525.1/)
10. *Brassica rapa*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000309985.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000309985.2/)
11. *Camelina sativa*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000633955.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633955.1/)
12. *Capsella rubella*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000375325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000375325.1/)
13. *Cardamine amara*
	- Assembly and annotation from Marek Šlenker (to be published)
14. *Cardamine glauca*
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (to be published)
	- Available at <https://github.com/mduchoslav/Genome_annotations_ultramafic_Brassicaceae>
15. *Cardamine hirsuta*
	- Gan, Xiangchao, Angela Hay, Michiel Kwantes, Georg Haberer, Asis Hallab, Raffaele Dello Ioio, Hugo Hofhuis, et al. “The Cardamine Hirsuta Genome Offers Insight into the Evolution of Morphological Diversity.” Nature Plants 2, no. 11 (October 31, 2016): 1–7. <https://doi.org/10.1038/nplants.2016.167>.
	- Data: <http://chi.mpipz.mpg.de/assembly.html>
16. *Cochlearia excelsa*
	- Bray, Sian M., Tuomas Hämälä, Min Zhou, Silvia Busoms, Sina Fischer, Stuart D. Desjardins, Terezie Mandáková, et al. “Kinetochore and Ionomic Adaptation to Whole-Genome Duplication in Cochlearia Shows Evolutionary Convergence in Three Autopolyploids.” Cell Reports 43, no. 8 (August 27, 2024). <https://doi.org/10.1016/j.celrep.2024.114576>.
	- Data: <https://doi.org/10.5061/dryad.ncjsxkt1s>
17. *Conringia planisiliqua*
	- Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig, Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria Willing, et al. “Improving and Correcting the Contiguity of Long-Read Genome Assemblies of Three Plant Species Using Optical Mapping and Chromosome Conformation Capture Data.” Genome Research 27, no. 5 (May 1, 2017): 778–86. <https://doi.org/10.1101/gr.213652.116>.
18. *Erysimum linariifolium*
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (to be published)
	- Available at <https://github.com/mduchoslav/Genome_annotations_ultramafic_Brassicaceae>
19. *Euclidium syriacum*
	- Jiao, Wen-Biao, Gonzalo Garcia Accinelli, Benjamin Hartwig, Christiane Kiefer, David Baker, Edouard Severing, Eva-Maria Willing, et al. “Improving and Correcting the Contiguity of Long-Read Genome Assemblies of Three Plant Species Using Optical Mapping and Chromosome Conformation Capture Data.” Genome Research 27, no. 5 (May 1, 2017): 778–86. <https://doi.org/10.1101/gr.213652.116>.
20. *Eutrema salsugineum*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000478725.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000478725.1/)
21. *Noccaea praecox*
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (to be published)
	- Available at <https://github.com/mduchoslav/Genome_annotations_ultramafic_Brassicaceae>
22. *Odontarrhena muralis*
	- Nezamivand-Chegini Mahnaz, Duchoslav Miloš, ..., Kolar Filip 2025 (to be published)
	- Available at <https://github.com/mduchoslav/Genome_annotations_ultramafic_Brassicaceae>
23. *Raphanus sativus*
	- Annotation from The NCBI Eukaryotic Genome Annotation Pipeline
	- [GCF_000801105.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000801105.1/)

### Removed from this run
24. *Rorippa islandica*
	- <https://phytozome-next.jgi.doe.gov/info/Rislandica_v1_1>

## Scripts and detailed description of methods

I use for documentation RMarkdown combining snippets of R code and Bash code. These Rmd files are then knitted to GitHub markdown files for better viewing.

0. SW installation (and description of versions used)
	- [GitHub md file](0_Installation_of_SW.md) | [original RMarkdown file](0_Installation_of_SW.Rmd)
1. Data preparation and [OrthoFinder](https://github.com/davidemms/OrthoFinder) run
	- [GitHub md file](1_Orthofinder_Brassicaceae.md) | [original RMarkdown file](1_Orthofinder_Brassicaceae.Rmd)
2. Supplementing orthologues in *A. thaliana* with other types of homologues
	- [GitHub md file](2_Supplementing_orthologues.md) | [original RMarkdown file](2_Supplementing_orthologues.Rmd)
3. Functional annotation of genes (from [InterProScan](https://interproscan-docs.readthedocs.io) and annotation of *A. thaliana* homologues)
	- [GitHub md file](3_Functional_annotations.md) | [original RMarkdown file](3_Functional_annotations.Rmd)
4. Statistics and plots (partially also in previous scripts)
	- [GitHub md file](4_Stats_and_plots.md) | [original RMarkdown file](4_Stats_and_plots.Rmd)
	
## Intermediate results

1. [Protein sequences of primary transcript](primary_transcripts/) - sequences of longest isoforms as input to OrthoFinder.
2. [OrthoFinder results](orthofinder_results/Results_brassicaceae_3/) - main files from OrthoFinder output.
3. [Length of protein sequences](protein_length_primary_transcripts/) - length of primary transcripts (in amino acid residues).
	- The same information is also in the final results.
4. [Supplemented orthologues](supplemented_orthologues/) - Orthologues between each species and *A. thaliana* from OrthoFinder supplemented with homologues using other methods.
	- Description of columns in the output tables is in [2_Supplementing_orthologues](2_Supplementing_orthologues.md#explanation-of-columns-in-output-table).
	- The same information is also in the final results.
	
Other intermediate results are not included as they are either not useful or they are too big (like the BLAST results or InterProScan results).

## Final results

The main results are tables (for each species except *A. thaliana*) with information for each gene including

- Protein lenght
- *A. thaliana* orthologues/homologues and details for their inference
- Gene ontology (GO) terms
- Other functional annotation from several sources

In case of several *A. thaliana* homologues assigned to one gene, annotation of all of them was compiled together.

There are several versions of tables (with different subsets of columns) in the [functional_annotation directory](functional_annotation/) and you can use the one which suits you:

1. [1_only_At_orthologues_and_GO](functional_annotation/1_only_At_orthologues_and_GO/)
	- Only the main columns: `Species_name` (gene IDs of the given species), `prot_length`, `Arabidopsis_thaliana`, `single_Arabidopsis_thaliana`, `homologue_type`, `GO_term_IDs` and `UniProt_Protein.names`
	- Format: plain tsv
2. [2_full_without_InterProScan_pathways_tsv_gz](functional_annotation/2_full_without_InterProScan_pathways_tsv_gz/)
	- All columns, just without `InterProScan_pathways` column, which contains a lot of data, making files too big.
	- Format: gzipped tsv
3. [3_InterProScan_pathways_tsv_gz](functional_annotation/3_InterProScan_pathways_tsv_gz/)
	- Only `Species_name` (gene IDs of the given species) and `InterProScan_pathways`
	- Format: gzipped tsv
4. [4_full_rds](functional_annotation/4_full_rds/)
	- All columns
	- Format: rds (can be imported to R using function `readRDS()`)
5. [5_full_with_support_for_my_annotations_rds](functional_annotation/5_full_with_support_for_my_annotations_rds/)
	- All columns plus gene model support columns (only for genomes annotated by me)
	- Format: rds (can be imported to R using function `readRDS()`)
	
The detailed description of the columns is in [3_Functional_annotations](3_Functional_annotations.md#explanation-of-columns-in-the-output-tables).
