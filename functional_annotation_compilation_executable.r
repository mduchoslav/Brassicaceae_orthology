### Extraction of gene/protein functional annotations from several sources

# Milos Duchoslav, 11/2025

# This script uses list of A. thaliana genes (AGI codes) and extracts various annotations for them.
# There can be several genes collapsed to one line. In such case, the script collapses also the annotations, if they differ (and can write there the ids for which the annotations apply).

# The script is heavily customised and hardly can be used universally.

## Usage of the script:
# Rscript --verbose "functional_annotation_compilation_executable.r" Species_name


### arguments from command line

# species
one.species <- commandArgs(trailingOnly = T)[1]


# chosen species for testing
# one.species <- "Cardamine_glauca"

# setwd("D:/!ecolgen/Brassicaceae_orthology/brassicaceae_3/")



### FUNCTIONS

# Sorry, this old function is described in Czech.

## Funkce k vytahování vlastností genů/proteinů z dalších databází
# Vytahuje informace pro skupiny genů (oddělené např. ";"). Pokud je informace pro celou skupinu stejná, je uvedena jen ta. Pokud je pro geny ze skupiny různá, uvede výčet a v závorce pro které geny.
# Funkce je míněna pro použití s funkcí sapply. Vstup do funkce (x) má být jeden textový řetězec (seznam genů uddělených oddělovačem). Pomocí funkce sapply lze použít na vektor takových textových řetězců.
# Parametry:
# x - textový řetězec s identifikátory genů oddělenými oddělovačem (viz "split")
# table - tabulka identifikátorů genů (resp. vektor) s přiřazenými vlastnostmi; pokud zde bude stejný gen vícekrát, nefunguje přiřazování genů k příslušným vlastnostem (gene.details)
# feature - vlastnosti genů (vektor, data.frame nebo matrix ve stejném pořadí jako "table")
# split - oddělovač použitý v "x" pro řetězení identifikátorů genů
# gene.details - TRUE/FALSE - Vyhodit v závorce, pro které geny platí jaká hodnota?
# collapse - znak, který bude vložen mezi hodnoty v případě více různých vlastností (feature) pro danou skupinu genů

## faster version of function

gene.properties <- function(x, table, feature, split = ", ", gene.details = T, collapse = ", ") {
	# rozdělit "x" podle "split" a vytáhnout k tomu vlastnosti (table - seznam identifikátorů; feature - vlastnosti přiřazené k identifikátorům)
	genes <- unlist(strsplit(x, split = split))
	row.numbers <- sapply(X = genes, FUN = grep, x = table, ignore.case = T)
	feature.v <- as.matrix(feature)[unlist(row.numbers), , drop = F]
	if(nrow(feature.v) == 1) {
		output <- feature.v
	} else if (nrow(feature.v) == 0) {
		output <- rep(NA, ncol(as.matrix(feature)))
	} else {
		output <- vector(mode = "character", length = ncol(as.matrix(feature)))
		for(j in 1:ncol(as.matrix(feature))) {
			feature.ind.v <- feature.v[, j]
			# Je "feature" u všech genů stejná?
			feature.identical <- all(feature.ind.v[1] == feature.ind.v)
			if(is.na(feature.identical)) {
				output[j] <- NA
			} else {
				if(feature.identical) {
					# v případě, že jsou všechny "feature" stejné, vyhoď tuto hodnotu
					output[j] <- feature.ind.v[1]
				} else {
					# v opačném případě načti ke každé hodnotě příslušné geny a ty vyhoď v závorce za každou hodnotou
					# (pokud gene.details = T)
					feature.f <- as.factor(feature.ind.v)
					if(gene.details) {
						feature.m <- matrix(nrow = nlevels(feature.f), ncol = 2)
						feature.m[, 1] <- levels(feature.f)
						for(i in 1:nlevels(feature.f)) {
							genes.level <- genes[feature.f == levels(feature.f)[i]]
							feature.m[i, 2] <- paste(genes.level, collapse = ";")
						}
						output[j] <- paste(paste0(feature.m[, 1], " (", feature.m[, 2], ")"), collapse = collapse)
					}
					else {
						# pokud gene.details = F, vyhoď hodnoty bez genů v závorce
						output[j] <- paste(levels(feature.f), collapse = collapse)
					}
				}
			}
		}
	}
	names(output) <- colnames(feature)
	return(output)
}


### Function for printing messages with time stamps

# Arguments:
# ... - Anything that will be pasted after the time stamp.
# sep - Separator of the arguments printed after the time stamp. Defaults to space.

t.cat <- function(..., sep = " ") {
	cat("[", format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "] ", paste(..., sep = sep), "\n", sep = "")
}







# Save console output to log file
dir.create(path = "functional_annotation/logs", recursive = T, showWarnings = F)
log.file <- file(paste0("functional_annotation/logs/", one.species, "_functional_annotation.log"), open = "wt")
sink(file = log.file, split = T)
sink(file = log.file, type = "message")


t.cat("Starting processing species", one.species)


### Reading common files

t.cat("Reading common files")

## Gene aliases
# File "gene_aliases_20241001.txt.gz"
# Downloaded on 2024-12-09 from
# https://v2.arabidopsis.org/download_files/Subscriber_Data_Releases/TAIR_Data_20240930/gene_aliases_20241001.txt.gz
# (web page https://v2.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSubscriber_Data_Releases%2FTAIR_Data_20240930)

aliases <- read.table(file = gzfile("functional_annotation/functional_annotation_data/Araport11/gene_aliases_20241001.txt.gz"), header = T, sep = "\t", quote = "", comment.char = "", fill = T)
# summary(aliases)
# dim(aliases)
# 
# head(aliases)
# tail(aliases)
# 
# # Are the gene IDs unique?
# table(aliases$locus_name)[table(aliases$locus_name)>1] # no
# aliases[aliases$locus_name == "AT1G01040", ]
# 
# ## there are some special characters that make problems
# # checking converting
# aa <- iconv(aliases$full_name, from = "UTF-8", to = "ASCII", sub = "__AAAA__")
# aa[grep("__AAAA__", aa)]
# aliases$full_name[grep("__AAAA__", aa)]

# remove special characters
aliases$full_name <- iconv(aliases$full_name, from = "UTF-8", to = "ASCII", sub = "")


## collapsing lines for the same genes
aliases.agg <- aggregate(x = aliases[, 2:3], by = list(aliases$locus_name), FUN = paste, collapse = ", ")
colnames(aliases.agg) <- colnames(aliases)
aliases.agg$full_name <- gsub(pattern = "NULL, |, NULL", replacement = "", x = aliases.agg$full_name)

# dim(aliases)
# dim(aliases.agg)


## Subcellular predictions
# File "Araport11-Subcellular_Predictions_version_2024_03_09.txt"
# Downloaded on 2024-12-09 from https://www.arabidopsis.org/download/list?dir=Genes%2FAraport11_genome_release

subcell <- read.table(file = "functional_annotation/functional_annotation_data/Araport11/Araport11-Subcellular_Predictions_version_2024_03_09.txt", header = F, sep = "\t")
# summary(subcell)
colnames(subcell) <- c("gene", "subcellular.prediction")

# # Are the gene IDs unique?
# table(subcell$gene)[table(subcell$gene)>1] # yes


## Functional descriptions from arabidopsis.org
# File "Araport11_functional_descriptions_20241001.txt.gz"
# Downloaded on 2024-12-09 from
# https://v2.arabidopsis.org/download_files/Subscriber_Data_Releases/TAIR_Data_20240930/Araport11_functional_descriptions_20241001.txt.gz 
# (web page https://v2.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSubscriber_Data_Releases%2FTAIR_Data_20240930)

fun.desc <- read.table(file = gzfile("functional_annotation/functional_annotation_data/Araport11/Araport11_functional_descriptions_20241001.txt.gz"), header = T, sep = "\t", quote = "", comment.char = "")
# summary(fun.desc)
# # Missing values are sometimes "NULL" and sometimes missing.
# colSums(fun.desc == "")
# colSums(fun.desc == "NULL")

# Making all missing values NULL
fun.desc[fun.desc == ""] <- "NULL"

# # Are the gene IDs unique?
# table(fun.desc$name)[table(fun.desc$name)>1] # yes

## GO annotations from arabidopsis.org
# Downloaded on 2025-11-11 from
# https://www.arabidopsis.org/download/file?path=Public_Data_Releases/TAIR_Data_20240930/ATH_GO_GOSLIM.txt.gz 
# (web page https://www.arabidopsis.org/download/list?dir=Public_Data_Releases%2FTAIR_Data_20240930)
go.ara <- read.table(file = unz(description = "functional_annotation/functional_annotation_data/Araport11/ATH_GO_GOSLIM.txt.gz", filename = "ATH_GO_GOSLIM.txt"), header = F, sep = "\t", quote = "", comment.char = "", skip = 4)

# summary(go.ara)
# dim(go.ara)

colnames(go.ara) <- c("locus", "acc_tair", "object_name", "rel_type", "GO_term", "GO_ID", "keyword_id_tair", 
											"aspect", "GOslim_term", "evid_code", "evid_descr", "evid_with", "ref", "annotator", "date")

# Column headers: explanation
# 1. locus name: standard AGI convention name. 
# 2. TAIR accession: the unique identifier for an object in the TAIR database -  the object type is the prefix, followed by a unique accession number(e.g. gene:12345). 
# 3. object name: the name of the object (gene, protein, locus) being annotated.
# 4. relationship type: the relationship between the annotated object and the GO term.
# 5. GO term: the actual string of letters corresponding to the GO ID.
# 6. GO ID: the unique identifier for a GO term. 
# 7. TAIR Keyword ID: the unique identifier for a keyword in the TAIR database.
# 8. Aspect: F=molecular function, C=cellular component, P=biological 13process.
# 9. GOslim term: high level GO term helps in functional categorization.
# 10. Evidence code: three letter code for evidence types (see: http://geneontology.org/page/guide-go-evidence-codes).
# 11. Evidence description: the analysis that was done to support the annotation.
# 12. Evidence with: supporting evidence for IGI, IPI, IC, IEA and ISS annotations.
# 13. Reference: Either a TAIR accession for a reference (reference table: reference_id) or reference from PubMed (e.g. PMID:1234). 
# 14. Annotator: TAIR, TIGR or a TAIR community member.
# 15. Date annotated: date the annotation was made.

# go.ara$locus
# # Are the gene IDs unique?
# table(go.ara$locus)[table(go.ara$locus)>1][1:60]

# collapsing lines for the same genes, choosing just GO_ID
go.ara.agg <- aggregate(x = go.ara[, c("locus", "GO_ID")], by = list(go.ara$locus), FUN = function(x) {paste(unique(x), collapse = ";")})[, -1]
# dim(go.ara.agg)


## Metabolic pathways
# File "aracyc_pathways.20230103"
# Downloaded on 2024-12-09 from
# https://plantcyc-ftp.storage.googleapis.com/pmn/Pathways/Data_dumps/PMN15.5_January2023/pathways/aracyc_pathways.20230103
# (web page https://plantcyc.org/)
plantcyc <- read.table("functional_annotation/functional_annotation_data/PlantCyc/aracyc_pathways.20230103", header = T, sep = "\t", quote = "", comment.char = "")
# colnames(plantcyc)
# head(plantcyc)
# # Are the gene IDs unique?
# table(plantcyc$Gene.id) # no

## collapsing lines for the same genes
plantcyc.agg <- aggregate(x = plantcyc, by = list(plantcyc$Gene.id), FUN = function(x) {paste(unique(x), collapse = ", ")})[, -1]


## Annotations from Uniprot
# File "uniprotkb_proteome_up000006548_2025_11_11.tsv.gz"
# Downloaded on 2025-11-11 from https://www.uniprot.org/uniprot/?query=proteome%3Aup000006548.
# I selected the desired columns for view and then downloaded that as a table.

uniprot.1 <- read.table(gzfile("functional_annotation/functional_annotation_data/Uniprot/uniprotkb_proteome_up000006548_2025_11_11.tsv.gz"), header = T, sep = "\t", quote = "", comment.char = "")
# summary(uniprot.1)
# head(uniprot.1)

# # Some proteins are annotated by several gene IDs. However, that should not be problem as the gene IDs are searched by grep.
# uniprot.1$Gene.Names..ordered.locus.[nchar(uniprot.1$Gene.Names..ordered.locus.) > 9][1:100]





t.cat("Starting reading specific files and processing them for", one.species)

## List of genes (AGI codes; there can be multiple genes per line separated by e.g. ";")

# data (gene numbers)
genes.pre0 <- read.table(file = paste0("supplemented_orthologues/", one.species, "__v__Arabidopsis_thaliana_supplemented.tsv"), sep = "\t", header = T)

genes <- genes.pre0[, c(1, 3)]
genes$ids <- genes$Arabidopsis_thaliana
# head(genes)

## 1a. Short version of gene IDs
short.ids.pre1 <- gsub(pattern = "\\.\\d", replacement = "", x = genes$ids)
short.ids.pre2 <- strsplit(short.ids.pre1, split = ", ")
short.ids.l <- lapply(X = short.ids.pre2, FUN = unique)
short.ids <- sapply(X = short.ids.l, FUN = paste, collapse = ", ")
short.ids[short.ids == "NA"] <- NA

genes$short.ids <- short.ids


## Gene aliases
t.cat("Extracting gene aliases for", one.species)

# extraction of gene properties by the function gene.properties
extr.aliases <- t(sapply(X = short.ids, FUN = gene.properties, table = aliases.agg$locus_name, 
												 feature = aliases.agg[, 2:3], gene.details = T, USE.NAMES = F))


## Subcellular predictions
t.cat("Extracting subcellular predictions for", one.species)

# extraction of gene properties by the function gene.properties
extr.subcell <- sapply(X = genes$ids, FUN = gene.properties, table = subcell[, 1], feature = subcell[, 2], USE.NAMES = F)


## Functional descriptions from arabidopsis.org
t.cat("Extracting functional descriptions from arabidopsis.org for", one.species)

# extraction of gene properties by the function gene.properties
extr.fun.desc <- t(sapply(X = genes$ids, FUN = gene.properties, table = fun.desc$name, feature = fun.desc[, 2:5], USE.NAMES = F))

## GO terms from arabidopsis.org
t.cat("Extracting GO terms from arabidopsis.org for", one.species)

# extraction of gene properties by the function gene.properties
extr.go.ara <- sapply(X = short.ids, FUN = gene.properties, table = go.ara.agg$locus, feature = go.ara.agg$GO_ID, gene.details = F, collapse = ";", USE.NAMES = F)

# each string split by ;, make the items unique and then again collapse them to one string
used.col <- sapply(X = extr.go.ara, FUN = function(x) {
	paste(unique(unlist(strsplit(x = x, split = ";"))), collapse = ";")
})
extr.go.ara <- used.col


## Metabolic pathways
t.cat("Extracting metabolic pathways for", one.species)

# extracting pathways
pathway <- sapply(X = short.ids, FUN = gene.properties, 
									table = plantcyc.agg$Gene.id, feature = plantcyc.agg$Pathway.name, gene.details = T, collapse = ", ", 
									USE.NAMES = F)


## Annotations from Uniprot
t.cat("Extracting annotations from Uniprot for", one.species)

# extraction of gene properties by the function gene.properties
# system.time(
uniprot.1.matched <- t(sapply(X = short.ids, FUN = gene.properties, 
															table = uniprot.1$Gene.Names..ordered.locus., 
															feature = uniprot.1, gene.details = F, collapse = ";", 
															USE.NAMES = F))
# )
# It can take an hour, would be good to somehow optimize.


## cleaning the matched Uniprot table

uniprot.1.matched.b <- uniprot.1.matched
# for each column
for(i in 1:ncol(uniprot.1.matched.b)) {
	used.col <- uniprot.1.matched.b[, i]
	# remove ; at line starts and ends
	used.col <- sub(pattern = "^;|;$", replacement = "", x = used.col)
	# remove double ;;
	used.col <- sub(pattern = ";;", replacement = ";", x = used.col)
	# save the column back
	uniprot.1.matched.b[, i] <- used.col
}

## make the Uniprot GO terms and similar columns unique for each protein of original species
# (removing duplicates in one cell)

# for each column
# some columns like "Cofactor" use ; in other way, it is better not to include them
for(i in c(17:27, 31:ncol(uniprot.1.matched.b))) {
	used.col <- uniprot.1.matched.b[, i]
	# remove ; at line starts and ends
	used.col <- sub(pattern = "^;|;$", replacement = "", x = used.col)
	# each string split by ; (can be surrounded by spaces), make the items unique and then again collapse them to one string
	used.col <- sapply(X = used.col, FUN = function(x) {
		paste(unique(unlist(strsplit(x = x, split = " *; *"))), collapse = ";")
	})
	
	uniprot.1.matched.b[, i] <- used.col
}

## remove unnecessary Uniprot columns
# colnames(uniprot.1.matched.b)
# columns to remove
cols.rem <- c("Gene.Names..ORF.", 
							"Organism",
							"Reactome",
							"PlantReactome",
							"AlphaFoldDB",
							"KEGG",
							"TAIR",
							"Araport"
)

uniprot.1.matched.c <- uniprot.1.matched.b[, !colnames(uniprot.1.matched.b) %in% cols.rem]


## Adding link to TAIR and Thalemine

# examples:
# https://bar.utoronto.ca/thalemine/gene:AT1G31690
# http://www.arabidopsis.org/servlets/TairObject?name=AT1G06550&type=locus

thalemine.link <- mapply(FUN = paste0, "https://bar.utoronto.ca/thalemine/gene:", short.ids.l, MoreArgs = list(collapse = "; "))
# remove in the rows with NA
thalemine.link[is.na(short.ids)] <- NA

tair.link <- mapply(FUN = paste0, "http://www.arabidopsis.org/servlets/TairObject?name=", short.ids.l, "&type=locus", MoreArgs = list(collapse = "; "))
tair.link[is.na(short.ids)] <- NA


## Completition of the table of annotations

annotations.1 <- cbind(genes.pre0, extr.aliases, extr.subcell, extr.fun.desc, extr.go.ara, pathway, uniprot.1.matched.c, thalemine.link, tair.link)
colnames(annotations.1)[colnames(annotations.1) == "extr.subcell"] <- "subcellular.prediction"
colnames(annotations.1)[colnames(annotations.1) == "extr.go.ara"] <- "GO"

# renaming columns (to reflect source)
colnames(annotations.1)[6:10] <- paste0("OrthoFinder_", colnames(annotations.1)[6:10])
colnames(annotations.1)[33:40] <- paste0("TAIR_", colnames(annotations.1)[33:40])
colnames(annotations.1)[41] <- paste0("PlantCyc_", colnames(annotations.1)[41])
colnames(annotations.1)[42:75] <- paste0("UniProt_", colnames(annotations.1)[42:75])



### InterproScan annotation of the original genes (not A. thaliana orthologues)

t.cat("Extracting InterproScan annotation for", one.species)

# read gzipped table
interpro <- read.table(file = gzfile(paste0("interproscan/", one.species, ".tsv.gz")), header = F, sep = "\t", quote = "", comment.char = "")

colnames(interpro) <- c("protein", "md5", "length", "analysis", "signature_accession", 
												"signature_description", "start", "stop", "score", "status", "run_date", 
												"interpro_accession", "interpro_description", "GO_term", "pathways")

# ## checks
# 
# # how many genes will be added?
# intepro.genes.uniq <- unique(interpro[, 1])
# length(intepro.genes.uniq)
# intepro.genes.new <- intepro.genes.uniq[! intepro.genes.uniq %in% annotations.1[, one.species]]
# length(intepro.genes.new)
# 
# interpro.new <- interpro[interpro[, 1] %in% intepro.genes.new, ]

# list of tools used by InterProScan
tools.interpro <- sort(unique(interpro$analysis))

# # How many rows and how many different values are there for each tool?
# for(i in tools.interpro) {
#   num.lines <- sum(interpro$analysis == i)
#   num.annot.prot <- length(unique(interpro$protein[interpro$analysis == i]))
#   nlev.tool <-  nlevels(as.factor(interpro$signature_accession[interpro$analysis == i]))
#   print(paste(i, "tool - records:", num.lines, "annotated proteins:", num.annot.prot, "levels:", nlev.tool))
# }
# 
# interpro[interpro$analysis == "AntiFam", ][1:10, ]

## making the interpro table wide (tools as separate columns)
# prepare the dataframe
wide.interpro <- data.frame(protein = unique(interpro[, 1]))
# loop through the tools
for(i in tools.interpro) {
	# extract rows for the tool
	extr.interpro <- interpro[interpro$analysis == i, ]
	# paste together description, accession and region
	extr.interpro.v <- paste0(extr.interpro$signature_description, " (", extr.interpro$signature_accession,
														"|", extr.interpro$start, "-", extr.interpro$stop, ")")
	# aggregate if there are several rows for one protein
	extr.interpro.agg <- aggregate(x = extr.interpro.v, by = list(extr.interpro$protein), FUN = paste, collapse = "; ")
	# extr.interpro.v[1:20]
	# extr.interpro.agg[1:20, ]
	
	# add as a new column
	wide.interpro[match(x = extr.interpro.agg[, 1], table = wide.interpro$protein), 
								which(tools.interpro == i) + 1] <- extr.interpro.agg[, 2]
	# rename the column
	colnames(wide.interpro)[which(tools.interpro == i) + 1] <- paste0("InterProScan_", i)
	
}


## add interpro accession and description

# extract rows with information
extr.interpro <- interpro[interpro$interpro_accession != "-", ]
# paste together description, accession and region
extr.interpro.v <- paste0(extr.interpro$interpro_description, " (", extr.interpro$interpro_accession,
													"|", extr.interpro$start, "-", extr.interpro$stop, ")")
# extr.interpro.v[1:30]
# aggregate if there are several rows for one protein
extr.interpro.agg <- aggregate(x = extr.interpro.v, by = list(extr.interpro$protein), FUN = paste, collapse = "; ")
# add as a new column
wide.interpro$InterProScan_interpro <- extr.interpro.agg[, 2][match(x = wide.interpro$protein, table = extr.interpro.agg[, 1])]


## add GO terms

# extract rows with information
extr.interpro <- interpro[interpro$GO_term != "-", ]
# remove origin of GO term in parentheses
extr.interpro.v <- gsub(pattern = "\\([[:alpha:]]*\\)", replacement = "", x = extr.interpro$GO_term)
# extr.interpro.v[1:30]
# aggregate if there are several rows for one protein (don't repeat the same GO terms)
extr.interpro.agg <- aggregate(x = extr.interpro.v, by = list(extr.interpro$protein), FUN = paste, collapse = "|")
# each string split by |, make the items unique and then again collapse them to one string
extr.interpro.agg$uniq <- sapply(X = extr.interpro.agg[, 2], FUN = function(x) {
	paste(unique(unlist(strsplit(x = x, split = "|", fixed = T))), collapse = ";")
})
# add as a new column
wide.interpro$InterProScan_GO_term <- extr.interpro.agg$uniq[match(x = wide.interpro$protein, table = extr.interpro.agg[, 1])]



## adding Pathways

# extract rows with information
extr.interpro <- interpro[interpro$pathways != "-", ]
extr.interpro.v <- extr.interpro$pathways
# extr.interpro.v[1:3]
# aggregate if there are several rows for one protein (don't repeat the same GO terms)
extr.interpro.agg <- aggregate(x = extr.interpro.v, by = list(extr.interpro$protein), FUN = paste, collapse = "|")
# extr.interpro.agg[1:3, ]
# each string split by |, make the items unique and then again collapse them to one string
extr.interpro.agg$uniq <- sapply(X = extr.interpro.agg[, 2], FUN = function(x) {
	paste(unique(unlist(strsplit(x = x, split = "|", fixed = T))), collapse = ";")
})
# extr.interpro.agg[1:3, ]
# add as a new column
wide.interpro$InterProScan_pathways <- extr.interpro.agg$uniq[match(x = wide.interpro$protein, table = extr.interpro.agg[, 1])]


## Integrating with annotations from A. thaliana

annotations.2 <- merge(x = annotations.1, y = wide.interpro, by.x = one.species, by.y = "protein", all = T)

# View(annotations.2[1:30, c("UniProt_Gene.Ontology.IDs", "InterProScan_GO_term")])

## Integrate A. thaliana (TAIR and Uniprot) and InterProScan GO terms (target species)

GO_term_IDs <- paste(annotations.2$TAIR_GO, annotations.2$"UniProt_Gene.Ontology.IDs", annotations.2$"InterProScan_GO_term", sep = ";")
# GO_term_IDs[1:10]

# remove NAs
GO_term_IDs.2 <- gsub(pattern = "NA", replacement = "", x = GO_term_IDs)
# GO_term_IDs.2[1:20]

# remove unnecessary ";"
GO_term_IDs.3 <- gsub(pattern = "^;*|;*$", replacement = "", x = GO_term_IDs.2)
# GO_term_IDs.3[1:10]

# each string split by "; ", make the items unique and then again collapse them to one string
GO_term_IDs.uniq <- sapply(X = GO_term_IDs.3, FUN = function(x) {
	paste(unique(unlist(strsplit(x = x, split = ";", fixed = T))), collapse = ";")
})

annotations.2$GO_term_IDs <- GO_term_IDs.uniq

# annotations.2$GO_term_IDs[1:10]


## reorder columns

# colnames(annotations.2)
annotations.3 <- annotations.2[, c(1:5, ncol(annotations.2), 6:(ncol(annotations.2)-1))]
# colnames(annotations.3)
# colnames(annotations.3) != "InterProScan_pathways"


### Exporting tables

t.cat("Exporting tables for", one.species)

# make folders
dir.create(path = "functional_annotation/1_only_At_orthologues_and_GO", recursive = T, showWarnings = F)
dir.create(path = "functional_annotation/2_full_without_InterProScan_pathways_tsv_gz", recursive = T, showWarnings = F)
dir.create(path = "functional_annotation/3_InterProScan_pathways_tsv_gz", recursive = T, showWarnings = F)
dir.create(path = "functional_annotation/4_full_rds", recursive = T, showWarnings = F)

## export full table without pathways (pathways are too big)

# gzipped table
write.table(x = annotations.3[, colnames(annotations.3) != "InterProScan_pathways"], 
						file = gzfile(paste0("functional_annotation/2_full_without_InterProScan_pathways_tsv_gz/", one.species, "_At_orthologues_and_functional_annotations.tsv.gz")), 
						sep = "\t", row.names = F, na = "", quote = F)

# export pathways
write.table(x = annotations.3[, c(one.species, "InterProScan_pathways")], 
						file = gzfile(paste0("functional_annotation/3_InterProScan_pathways_tsv_gz/", one.species, "_InterProScan_pathways.tsv.gz")), 
						sep = "\t", row.names = F, na = "", quote = F)

# write RDS
saveRDS(object = annotations.3, 
				file = paste0("functional_annotation/4_full_rds/", one.species, "_At_orthologues_and_functional_annotations_full.rds"))


## export main columns only
write.table(x = annotations.3[, c(1:6, 45)], 
						file = paste0("functional_annotation/1_only_At_orthologues_and_GO/", one.species, "_At_orthologues_and_GO.tsv"),
						sep = "\t", row.names = F, na = "", quote = F)

t.cat(one.species, "finished.")

# Stop saving output to log file
sink(type = "message")
sink()
close(log.file)