#!/usr/bin/env Rscript

## PACKAGES ----
required_pkgs <- c("stringr", "microseq", "tidyr", "janitor", "optparse")

missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  message("\n\t", "Installing missing R packages...")
  install.packages(missing_pkgs)
}

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))

## OPTIONS ----
option_list <- list(
  make_option(c("-f", "--fasta"), 
              type = "character", 
              help = "Input sequence file (FASTA)"),
  
  make_option(c("-i", "--in_dir"), 
              type = "character", 
              help = "Directory containing VIGOR results")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if in_dir is provided
if (is.null(opt$fasta)) {
  print_help(opt_parser)
  stop("Input sequence file (fasta) is required.")
}

# Check if in_dir is exists
if (!file.exists(opt$fasta)) {
  print_help(opt_parser)
  stop("Input sequence file (fasta) does not exist.")
}

# Check if in_dir is provided
if (is.null(opt$in_dir)) {
  print_help(opt_parser)
  stop("Input directory (in_dir) is required.")
}

# Check if in_dir is exists
if (!dir.exists(opt$in_dir)) {
  print_help(opt_parser)
  stop("Input directory (in_dir) is does not exist.")
}

## IN_DIR & IN_FASTA ----
in_fasta <- opt$fasta
in_dir <- opt$in_dir

# in_dir <- "~/projects/rota-annotate/data/TEST_DATA/2_analysis_input/1/vigor"
# in_fasta <- "~/projects/rota-annotate/data/TEST_DATA/2_analysis_input/1/Assembly_1-3101566450_S2.fas"

## FUNCTIONS ----
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

concat_fas <- function(file_list, add_file_name, full_name){
  fas <- data.frame()
  for (i in 1:length(file_list)) {
    if (file.size(file_list[i]) == 0) {
      cat("concat_fas: following file is empty\n",
          paste0(" file_list[", i, "]: "), file_list[i])
      next(i)
    } else {
      df <- readFasta(file_list[i])
      if(isTRUE(add_file_name)) {
        if (isTRUE(full_name)) {
          df$file_name <- file_list[i]
        }
        if(isFALSE(full_name)) {
          df$file_name <- basename(file_list[i])
        }
      }
      fas <- rbind(fas, df)
    }
  }
  return(fas)
}

build_fas <- function(Header, Sequence) {
  fas <- data.frame(Header, Sequence, stringsAsFactors = F)
  return(fas)
}

concat_tables <- function(file_list, header, add_file_name, full_name) {
  df <- data.frame()
  for (i in 1:length(file_list)) {
    sub_df <- read.table(file = file_list[i], header = header, sep = '\t', comment.char = "#")
    if(isTRUE(add_file_name)) {
      if (isTRUE(full_name)) {
        sub_df$file_name <- file_list[i]
      }
      if(isFALSE(full_name)) {
        sub_df$file_name <- basename(file_list[i])
      }
    }
    df <- rbind(df, sub_df)
  }
  return(df)
}

split_header <- function(string, pattern) {
  pattern_1 <- paste0(pattern, "=.+$")
  pattern_2 <- paste0(pattern, '=\\"')
  val <- str_extract(string = string, pattern = pattern_1) %>% 
    str_remove(pattern = pattern_2) %>% 
    str_remove(pattern = '\\".+$') %>% 
    str_remove(pattern = '\\"$')
  return(val)
}

rename_col <- function(df, old_name, new_name) {
  if (isTRUE(old_name %in% colnames(df))) {
    names(df)[names(df) == old_name] <- new_name
    return(df)
  } else {
    cat("rename_col: old_name='", old_name, "', not found in colnames.", sep = "")
    return(df)
  }
}

read_RPT <- function(in_rpt) {
  rpt_out <- readLines(in_rpt)
  rpt_out <- rpt_out[nzchar(trimws(rpt_out))]
  
  ## split at lines that begin with "Sequence"
  idx <- grep(pattern = "^Sequence", x = rpt_out, value = FALSE)
  rpt_out <- splitAt(x = rpt_out, pos = idx)[-1]
  return(rpt_out)
}

parse_RPT <- function(rpt_in) {
  rpt_out <- list()
  for (i in 1:length(rpt_in)) {
    idx <- grep(pattern = "^(gene_id|Getting)", 
                x = rpt_in[[i]], 
                value = FALSE)
    sub_rpt <- splitAt(x = rpt_in[[i]], pos = idx)
    
    sub_rpt <- sub_rpt[1:2] 
    SEQUENCE <- str_split_fixed(string = sub_rpt[[1]], pattern = "\\s+", n = 10) %>% as.data.frame()
    SEQUENCE[1, ] <- SEQUENCE[1, ] %>% str_remove(pattern = "\\%") %>% str_trim()
    SEQUENCE <- row_to_names(dat = SEQUENCE, row_number = 1)
    sub_rpt[[1]] <- SEQUENCE
    names(sub_rpt)[1] <- "SEQUENCE"
    
    ## Gene
    GENE <- sub_rpt[[2]]
    idx <- str_which(string = GENE, pattern = "^\\s+", negate = TRUE)
    GENE <- GENE[idx]
    GENE <- str_split_fixed(string = GENE, pattern = "\\s+", n = 8)[,1:7] %>% as.data.frame()
    GENE[1, ] <- GENE[1, ] %>% str_remove(pattern = "\\%") %>% str_trim()
    GENE <- row_to_names(dat = GENE, row_number = 1)
    sub_rpt[[2]] <- GENE
    
    names(sub_rpt)[2] <- "GENE"
    rpt_out[[i]] <- sub_rpt
    names(rpt_out)[i] <- sub_rpt$SEQUENCE$Sequence[1]
  }
  return(rpt_out)
}


## FILE_LIST ----
file_list <- list.files(path = in_dir, 
                        full.names = TRUE, recursive = TRUE)

## ALN ----
# in_aln <- grep(pattern = "\\.aln$", x = file_list, value = TRUE)
# aln <- readLines(in_aln)
# idx <- str_which(string = aln, pattern = "^C4 Alignment:")
# aln <- splitAt(x = aln, pos = idx)

## RPT ----
in_rpt <- grep(pattern = "\\.rpt$", 
               x = file_list, value = TRUE)

# Check if the file exists
if (!file.exists(in_rpt)) {
  stop("File does not exist.")
}

rpt <- read_RPT(in_rpt = in_rpt)
rpt <- parse_RPT(rpt_in = rpt)

## RPT seq & gene ----
rpt_seq <- rpt %>% 
  lapply("[[", 'SEQUENCE') %>% 
  bind_rows()

rpt_seq$Ref_DB <- basename(rpt_seq$Ref_DB)

rpt_gene <- rpt %>% 
  lapply("[[", 'GENE') %>% 
  bind_rows()


colnames(rpt_gene) <- c("ORF_ID", "Percent_Identity", "Percent_Similarity",
                        "Percent_Coverage", "T5", "Gap", "T3")

## CDS ----
in_cds <- grep(pattern = "\\.cds$", 
               x = file_list, value = TRUE)

cds <- concat_fas(file_list = in_cds, 
                  add_file_name = TRUE, full_name = FALSE)

## CDS Seq_ID ----
cds$Seq_ID <- str_extract(string = cds$Header, 
                          pattern = "^[[:alnum:]_]+")

## CDS ORF_ID ----
cds$ORF_ID <- str_extract(string = cds$Header, 
                          pattern = "^\\S*")

## GFF_CDS ----
# when location is on reverse strand, 
# there is an issue with cds generated by vigor4
# cds and pep now retrieved by gffread ussing gff3 file
in_gff_cds <- grep(pattern = "GFF3_CDS\\.fas$",
                   x = file_list, value = TRUE)

gff_cds <- concat_fas(file_list = in_gff_cds, 
                      add_file_name = FALSE, full_name = FALSE)

colnames(gff_cds) <- c("ORF_ID", "ORF_NT")

col_names <- union(names(cds), names(gff_cds))

cds <- merge(x = cds, 
             y = gff_cds, 
             by = "ORF_ID", sort = FALSE)[, col_names]

cds$Sequence <- NULL

## PEP ----
in_pep <- grep(pattern = "\\.pep$", 
               x = file_list, value = TRUE)

pep <- concat_fas(file_list = in_pep, 
                  add_file_name = F, full_name = FALSE)

pep$Header <- str_extract(string = pep$Header, 
                          pattern = "^\\S*")

pep <- rename_col(df = pep,
                  old_name = "Header",
                  new_name = "ORF_ID")

## GFF_PEP ----
in_gff_pep <- grep(pattern = "GFF3_PEP\\.fas$",
                   x = file_list, value = TRUE)

gff_pep <- concat_fas(file_list = in_gff_pep, 
                      add_file_name = FALSE, full_name = FALSE)

colnames(gff_pep) <- c("ORF_ID", "ORF_AA")

col_names <- union(names(pep), names(gff_pep))

pep <- merge(x = pep, 
             y = gff_pep, 
             by = "ORF_ID", sort = FALSE)[, col_names]

pep$Sequence <- NULL

## MERGE CDS & PEP ----
cds <- full_join(x = cds, y = pep, by = "ORF_ID")

## CDS Location ----
cds$Location <- str_extract(string = cds$Header, 
                            pattern = "location=\\S*") %>% 
  str_remove(pattern = "location=")

cds$Location_Start <- str_extract(string = cds$Location, pattern = '^<*[0-9]+') %>% 
  str_extract(pattern = "[0-9]+$") %>% 
  as.numeric()

cds$Location_Stop <- str_extract(string = cds$Location,pattern = '>*[0-9]+$') %>% 
  str_extract(pattern = "[0-9]+$") %>% 
  as.numeric()

## CDS Codon_Start ----
cds$Codon_Start <- str_extract(string = cds$Header, pattern = "codon_start=\\S*") %>% 
  str_remove(pattern = "codon_start=")

## CDS Gene ----
cds$Gene <- split_header(string = cds$Header, 
                         pattern = "gene")

## CDS Product ----
cds$Product <- split_header(string = cds$Header, 
                            pattern = "product")

## CDS Ref_DB ----
cds$Ref_DB <- split_header(string = cds$Header, 
                           pattern = "ref_db")

## CDS Ref_ID ----
cds$Ref_ID <- split_header(string = cds$Header, 
                           pattern = "ref_id")

## CDS Length_ORF_NT ----
cds$Length_ORF_NT <- str_length(string = cds$ORF_NT)

## CDS Length_ORF_AA ----
cds$Length_ORF_AA <- str_length(string = cds$ORF_AA)

## CDS Merge rpt_gene and CDS ----
by <- intersect(x = colnames(rpt_gene), 
                y = colnames(cds))

rpt_gene <- full_join(x = rpt_gene, 
                      y = cds, 
                      by = by)

rpt_gene$seq.ID <- str_remove(string = rpt_gene$ORF_ID, 
                              pattern = "\\.[:alnum:]+$")

## FASTA ----
fasta <- readFasta(in.file = in_fasta)
colnames(fasta) <- c("seq.ID", "Sequence_forward")

## FASTA Merge rpt_gene and fasta ----
rpt_gene <- merge(x = rpt_gene, 
                  y = fasta, 
                  by = "seq.ID", all.x = TRUE, sort = FALSE)

idx <- which(rpt_gene$Location_Start > rpt_gene$Location_Stop)

rpt_gene$Sequence_forward[idx] <- reverseComplement(rpt_gene$Sequence_forward[idx])

## ORF Complete ----
rpt_gene$ORF <- NA
idx <- str_which(string = rpt_gene$Location, 
                          pattern = "^[0-9]+\\.\\.[0-9]+$")

rpt_gene$ORF[idx] <- 'complete'

## ORF Partial ----
idx <- str_which(string = rpt_gene$Location, 
                        pattern = "^[0-9]+\\.\\.[0-9]+$", 
                        negate = T)

rpt_gene$ORF[idx] <- 'partial'


## NOTES Gaps/frame shifts ----

rpt_gene$Notes <- ''

idx <- str_which(string = rpt_gene$Location, 
                 pattern = ",")

if (length(idx) > 0) {
  rpt_gene$Notes[idx] <- "Contains gap or frasmeshift"
}

## NOTES Multiple ORF ----

idx <- str_which(string = rpt_gene$ORF_ID, 
                 pattern = "\\.[0-9]+[a-z]+$")

if (length(idx) > 0) {
  rpt_gene$Notes[idx] <- paste0(rpt_gene$Notes[idx], "; Contains multiple ORFs")
}

rpt_gene$Notes <- str_remove(string = rpt_gene$Notes, 
                             pattern = "^; ")

rpt_gene[rpt_gene == ""] <- NA

## SUMMARY rpt_seq ----
rpt_seq <- rename_col(df = rpt_seq, 
                      old_name = "Sequence", new_name = "Seq_ID")

rpt_seq <- rename_col(df = rpt_seq, 
                      old_name = "Length", new_name = "Length_Seq_ID")

out_rpt_seq <- paste0(in_rpt, "_seq")

write.table(x = rpt_seq, 
            file = out_rpt_seq, 
            sep = '\t', na = 'NA', 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## SUMMARY rpt_gene ----

col_order <-  c("Seq_ID", "ORF_ID", "Notes", "ORF", "Gene", "Product", "Location", 
                "Codon_Start", "Length_ORF_NT", "Length_ORF_AA", "Ref_DB", 
                "Ref_ID", "Percent_Identity", "Percent_Similarity", 
                "Percent_Coverage", "T5", "Gap", "T3", "file_name", 
                "Sequence_forward", "ORF_NT", "ORF_AA")

rpt_gene <- rpt_gene[col_order]

out_rpt_gene <- paste0(in_rpt, "_gene")

write.table(x = rpt_gene, 
            file = out_rpt_gene, 
            sep = '\t', na = 'NA', 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
