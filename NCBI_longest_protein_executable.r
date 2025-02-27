### Extraction of the longest protein per gene from genomes annotated by NCBI Eukaryotic Genome Annotation Pipeline

# The script extracts information from GFF annotation file which protein belongs to which transcript and which gene. It uses this information to extract only the longest protein per gene from the fasta file of all proteins of an organism.

# Output:
# Fasta file with the longest protein per gene.
# Table listing protein IDs, their mRNA parents and their gene parents.

# Usage:
# The script can be run like this:
# Rscript --verbose "NCBI_longest_protein_executable.r" input_fasta_file input_gff_file


# Miloš Duchoslav (duchmil@gmail.com), 2022-2023




## Three functions copied from the package "seqinr"
# The functions are needed for independency on installed packages.
read.fasta <- function (file = system.file("sequences/ct.fasta.gz", package = "seqinr"), 
                       seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE, 
                       set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, 
                       strip.desc = FALSE, whole.header = FALSE, bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong, 
                       endian = .Platform$endian, apply.mask = TRUE) {
  seqtype <- match.arg(seqtype)
  if (!bfa) {
    lines <- readLines(file)
    if (legacy.mode) {
      comments <- grep("^;", lines)
      if (length(comments) > 0) 
        lines <- lines[-comments]
    }
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseq <- length(ind)
    if (nseq == 0) {
      stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                         collapse = ""))
    if (seqonly) 
      return(sequences)
    if (!whole.header) {
      nomseq <- lapply(seq_len(nseq), function(i) {
        firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
        substr(firstword, 2, nchar(firstword))
      })
    }
    else {
      nomseq <- lapply(seq_len(nseq), function(i) {
        substr(lines[ind[i]], 2, nchar(lines[ind[i]]))
      })
    }
    if (seqtype == "DNA") {
      if (forceDNAtolower) {
        sequences <- as.list(tolower(sequences))
      }
    }
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        Annot <- lines[ind[i]]
        if (strip.desc) 
          Annot <- substr(Annot, 2L, nchar(Annot))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = switch(seqtype, AA = "SeqFastaAA", 
                                                                         DNA = "SeqFastadna"))
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
  if (bfa) {
    if (seqtype != "DNA") 
      stop("binary fasta file available for DNA sequences only")
    mycon <- file(file, open = "rb")
    r2s <- words(4)
    readOneBFARecord <- function(con, sizeof.longlong, endian, 
                                 apply.mask) {
      len <- readBin(con, n = 1, what = "int", endian = endian)
      if (length(len) == 0) 
        return(NULL)
      name <- readBin(con, n = 1, what = "character", endian = endian)
      ori_len <- readBin(con, n = 1, what = "int", endian = endian)
      len <- readBin(con, n = 1, what = "int", endian = endian)
      seq <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                     size = 1, endian = endian)
      mask <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                      size = 1, endian = endian)
      if (endian == "little") {
        neword <- sizeof.longlong:1 + rep(seq(0, (len - 
                                                    1) * sizeof.longlong, by = sizeof.longlong), 
                                          each = sizeof.longlong)
        seq <- seq[neword]
        mask <- mask[neword]
      }
      seq4 <- c2s(r2s[as.integer(seq) + 1])
      seq4 <- substr(seq4, 1, ori_len)
      if (apply.mask) {
        mask4 <- c2s(r2s[as.integer(mask) + 1])
        mask4 <- substr(mask4, 1, ori_len)
        npos <- gregexpr("a", mask4, fixed = TRUE)[[1]]
        for (i in npos) substr(seq4, i, i + 1) <- "n"
      }
      return(list(seq = seq4, name = name))
    }
    sequences <- vector(mode = "list")
    nomseq <- vector(mode = "list")
    i <- 1
    repeat {
      res <- readOneBFARecord(mycon, sizeof.longlong, endian, 
                              apply.mask)
      if (is.null(res)) 
        break
      sequences[[i]] <- res$seq
      nomseq[[i]] <- res$name
      i <- i + 1
    }
    close(mycon)
    nseq <- length(sequences)
    if (seqonly) 
      return(sequences)
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        if (!strip.desc) 
          Annot <- c2s(c(">", nomseq[[i]]))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = "SeqFastadna")
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
}


write.fasta <- function (sequences, names, file.out, open = "w", nbchar = 60, 
                         as.string = FALSE) {
  outfile <- file(description = file.out, open = open)
  write.oneseq <- function(sequence, name, nbchar, as.string) {
    writeLines(paste(">", name, sep = ""), outfile)
    if (as.string) 
      sequence <- s2c(sequence)
    l <- length(sequence)
    q <- floor(l/nbchar)
    r <- l - nbchar * q
    if (q > 0) {
      sapply(seq_len(q), function(x) writeLines(c2s(sequence[(nbchar * 
                                                                (x - 1) + 1):(nbchar * x)]), outfile))
    }
    if (r > 0) {
      writeLines(c2s(sequence[(nbchar * q + 1):l]), outfile)
    }
  }
  if (!is.list(sequences)) {
    write.oneseq(sequence = sequences, name = names, nbchar = nbchar, 
                 as.string = as.string)
  }
  else {
    n.seq <- length(sequences)
    sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]), 
                                                    name = names[x], nbchar = nbchar, as.string = as.string))
  }
  close(outfile)
}

c2s <- function (chars = c("m", "e", "r", "g", "e", "d")) {
  return(paste(chars, collapse = ""))
}


## NCBI translator - translation table from protein ID to mRNA ID and to gene ID

# arguments from command line
arg.fasta <- commandArgs(trailingOnly = T)[1]
arg.gff <- commandArgs(trailingOnly = T)[2]


# arg.fasta = "D:/!ecolgen/resources/orthofinder/arenosa_lyrata_thaliana/A_lyrata_NCBI101.fasta"
# arg.fasta = "C_rubella.fasta"
# arg.fasta = "B_oleracea.fasta"
# arg.gff = "B_oleracea.gff"
# setwd("D:/!ecolgen/resources/orthofinder/brassicaceae_1/")


# read gff file
# gff <- read.table(file = "D:/!ecolgen/resources/ncbi_lyrata_reference_101/GCF_000004255.2_v.1.0_genomic.gff",
#                   header = F, sep = "\t", comment.char = "#", quote = ""
#                   # , nrows = 1500
# )
# gff <- read.table(file = "C_rubella.gff",
#                   header = F, sep = "\t", comment.char = "#", quote = ""
#                   # , nrows = 1500
# )
gff <- read.table(file = arg.gff,
                  header = F, sep = "\t", comment.char = "#", quote = "" 
                  # , nrows = 1500
)
cat("GFF", arg.gff, "with ", nrow(gff), "rows and", ncol(gff), "columns was read.\n\n")
# extract rows with CDS
gff.cds <- gff[gff$V3 == "CDS", ]
# extract ID
gff.cds$id <- gsub(pattern = "(^ID=cds-)|(;.*$)", replacement = "", x = gff.cds$V9)
# extract parent link
gff.cds$parent <- gsub(pattern = "(.*Parent=rna-)|(;.*$)", replacement = "", x = gff.cds$V9)

# number of lines (CDS) that does not have as parent RNA (e.g. chloroplast genes)
cds.without.rna.parent <- length(grep(pattern = "Parent=rna", x = gff.cds$V9, invert = T))
cat("There are", cds.without.rna.parent, "CDS that does not have RNA as the parent (counted as CDS lines in GFF).\nThese might be e.g. chloroplast or mitochondrial genes. Check that if the number is higher than 0!\n")


# Special case are chloroplast genes because they does not have mRNA (at least in A. lyrata).
# It seems that mitochondrial genes have mRNA annotated, at least in B. oleracea.
# View(gff.cds[grep(pattern = "^YP_", x = gff.cds$id), ])
gff.cds$parent[grepl(pattern = "^YP_", x = gff.cds$id) & grepl(pattern = "Parent=gene-", x = gff.cds$V9)] <- 
  gsub(pattern = "(.*Parent=gene-)|(;.*$)", replacement = "", x = gff.cds$V9[grep(pattern = "^YP_", x = gff.cds$id)])
cat("Parents in", 
    length(gff.cds$parent[grepl(pattern = "^YP_", x = gff.cds$id) & grepl(pattern = "Parent=gene-", x = gff.cds$V9)]),
    "chloroplast genes were corrected.\n\n")

# take only rows with unique id - parent combination
gff.cds.uniq <- gff.cds[!duplicated(gff.cds[, 10:11]), ]
cat("There are", nrow(gff.cds.uniq), "CDS annotated in the GFF.\n")
# nrow(gff.cds.uniq)

# extract rows with mRNA
gff.mrna <- gff[gff$V3 == "mRNA", ]
cat("There are", nrow(gff.mrna), "mRNA annotated in the GFF.\n")

# extract ID
gff.mrna$id <- gsub(pattern = "(^ID=rna-)|(;.*$)", replacement = "", x = gff.mrna$V9)
# extract parent link
gff.mrna$parent <- gsub(pattern = "(.*Parent=gene-)|(;.*$)", replacement = "", x = gff.mrna$V9)
# take only rows with unique id - parent combination
gff.mrna.uniq <- gff.mrna[!duplicated(gff.mrna[, 10:11]), ]
cat("There are", nrow(gff.mrna.uniq), "unique mRNA and their parent gene combinations in the GFF.\n\n")


# make the translation table
ncbi.translator <- data.frame(protein = gff.cds.uniq$id, 
                              mRNA = gff.cds.uniq$parent, 
                              gene = gff.mrna.uniq$parent[match(gff.cds.uniq$parent, gff.mrna.uniq$id)])
# for chloroplast genes, transfer the info from mRNA to gene column
ncbi.translator$gene[grepl(pattern = "^YP_", x = ncbi.translator$protein) & is.na(ncbi.translator$gene)] <- 
  ncbi.translator$mRNA[grepl(pattern = "^YP_", x = ncbi.translator$protein) & is.na(ncbi.translator$gene)]
# ncbi.translator$mRNA[grepl(pattern = "^YP_", x = ncbi.translator$protein) & ncbi.translator$mRNA == ncbi.translator$gene] <- NA

# name of the input fasta file
name.fasta.in <- sub(pattern = "^.*/", replacement = "", x = arg.fasta)
# adding of suffix to the output table file
name.table.out <- sub(pattern = "\\..*$", replacement = "_gene_mRNA_protein_translation_table.tsv", x = name.fasta.in)
# creation of directory for translation tables if it does not exist
if(!dir.exists("gene_mRNA_protein_translation_tables")) dir.create("gene_mRNA_protein_translation_tables")
# exporting the translation table
write.table(x = ncbi.translator,
            file = paste0("gene_mRNA_protein_translation_tables/", name.table.out),
            sep = "\t", quote = F, row.names = F)
cat("Gene - mRNA - protein translation table was saved under name:\n", paste0("gene_mRNA_protein_translation_tables/", name.table.out), "\n\n")


## reading the fasta file with all protein sequences
# library(seqinr)
fasta.ncbi <- read.fasta(file = arg.fasta, as.string = T)
# fasta.ncbi <- read.fasta(file = "D:/!ecolgen/resources/orthofinder/arenosa_lyrata_thaliana/A_lyrata_NCBI101.fasta", as.string = T)
cat("Fasta file", arg.fasta, "with", length(fasta.ncbi), "sequences was read.\n")
# names of sequences in fasta file (protein names)
fasta.names <- attr(x = fasta.ncbi, which = "name")
# names of corresponding genes
gene.names <- ncbi.translator$gene[match(x = fasta.names, table = ncbi.translator$protein)]
cat("There are", sum(duplicated(gene.names)), "alternative transcripts that should be removed.\n\n")

# length of sequences
len.seq <- sapply(X = fasta.ncbi, FUN = nchar)

# make a data.frame and order it according to gene names and length
fasta.df <- cbind.data.frame(fasta.names, gene.names, len.seq)[order(gene.names, len.seq, decreasing = T), ]
# for the rows with the same gene name, mark only the first one as TRUE (inverted duplicated) - because of the sorting it should be the longest protein
fasta.df$longest.transcript <- !duplicated(fasta.df$gene.names)
# filter for longest transcripts
fasta.df.longest <- fasta.df[fasta.df$longest.transcript, ]
# sanity checks
cat("Sanity checks:\n")
cat("Number of sequences to export:", nrow(fasta.df.longest), "\n")
cat("Unique gene names that are parents of the proteins in fasta file:", length(unique(gene.names)), "\n")
cat("Total number of sequences in fasta file (", length(fasta.ncbi), ") without number of alternative transcripts (",
    sum(duplicated(gene.names)), ") is: ", length(fasta.ncbi)-sum(duplicated(gene.names)), "\n\n", sep = "")


# filter the fastas
fasta.longest <- fasta.ncbi[names(fasta.ncbi) %in% fasta.df.longest$fasta.names]
# names of sequences in fasta file (protein names)
fasta.longest.names <- attr(x = fasta.longest, which = "name")

# adding of suffix to the output fasta file
name.fasta.out <- sub(pattern = "\\..*$", replacement = ".fasta", x = name.fasta.in)
# write fasta
# write.fasta(sequences = fasta.longest, names = fasta.longest.names, file.out = "D:/!ecolgen/resources/orthofinder/arenosa_lyrata_thaliana/A_lyrata_NCBI101_longest_transcript_2.fasta")
# creation of directory for exported fasta if it does not exist
if(!dir.exists("longest_transcripts")) dir.create("longest_transcripts")
write.fasta(sequences = fasta.longest, names = fasta.longest.names, file.out = paste0("longest_transcripts/", name.fasta.out))
cat("Fasta file with longest transcripts per gene (primary transcripts) was saved under name:\n", 
    paste0("longest_transcripts/", name.fasta.out), "\n")
