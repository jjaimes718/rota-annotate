#!/usr/bin/Rscript

## PACKAGES ----
required_pkgs <- c("stringr", "microseq", "openxlsx", "tools", "gtools", "optparse", "plyr")

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

# OPTIONS ----
option_list <- list(
  make_option(c("-i", "--in_dir"), 
              type = "character", 
              help = "Directory containing assembly results")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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

## IN_DIR ----
in_dir <- opt$in_dir

#in_dir <- "~/projects/rota-annotate/data/TEST_DATA"
in_dir <- normalizePath(path = in_dir, mustWork = TRUE)

project_name <- basename(in_dir)

message('RUNNING: 4_generate_analysis_summary.R', '\n\n\t', 
        'IN_DIR: ', in_dir)

## FUNCTIONS ----
build_fas <- function(Header,Sequence) {
  fas <- data.frame(Header, Sequence, stringsAsFactors = F)
  return(fas)
}

merge_dfs <- function(x,y,by,all) {
  df <- merge(x = x, y = y, by = by, all = all, sort = F)[, union(names(x), names(y))]
  return(df)
}

get_top_x <- function(df, x, grp_by_col) {
  top_x <- df %>% group_by(df[grp_by_col]) %>% slice(1:x)
  return(top_x)
}

get_sacc <- function(df, sacc_col) {
  idx <- which(!duplicated(df[sacc_col]))
  sacc <- df[idx,]
  sacc <- sacc[sacc_col]
  colnames(sacc) <- 'sacc'
  return(sacc)
}

concat_tables <- function(file_list, header) {
  size <- file.size(file_list)
  idx <- which(size > 0)
  file_list <- file_list[idx]
  df <- do.call("rbind", lapply(file_list, function(fn) 
    data.frame(read.table(file = fn, header = header, sep = '\t', fill = T))
  ))
  return(df)
}

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

get_idx_dup <- function(in_vector){
  idx <- which(duplicated(in_vector) | duplicated(in_vector, fromLast = TRUE))
  return(idx)
}

get_gtypes <- function(df, gt_colname ,segment) {
  test_seg <- segment %in% df$vigor_Gene
  
  if (isFALSE(test_seg)) {
    genotype <- order$gt[order$vigor_Gene == segment]
    n = 0
    final_gtype <- paste0(genotype,'(', n, ')')
  } else {
    idx <- which(df$vigor_Gene == segment)
    sub_seg <- df[idx,]
    if (any(is.na(sub_seg[[gt_colname]]))) {
      idx <- which(is.na(sub_seg[[gt_colname]]))
      
    }
    freq <- count(sub_seg, 
                  vars = c(gt_colname, "orf.2"))
    freq <- freq[mixedorder(freq[[gt_colname]]),]
    
    freq_gt <- count(df = freq[c(gt_colname)], 
                     vars = gt_colname)
    
    gt <- unique(freq_gt[[gt_colname]])
    
    for (y in 1:length(gt)){
      if(any((freq$orf.2 == 'complete') & (freq[[gt_colname]] == gt[y]))) {
        idx <- which(freq_gt[[gt_colname]] == gt[y])
        freq_gt[[gt_colname]][[idx]] <- paste0(freq_gt[[gt_colname]][[idx]], '*')
      }
    }
    
    final_gtype <- paste0(freq_gt[[gt_colname]], '(',freq_gt$freq, ')')
    final_gtype <- paste0(final_gtype, collapse = '')
  }
  final_gtype <- str_remove_all(string = final_gtype, pattern = "\\(1\\)")
  return(final_gtype)
}

get_gtype_simple <- function(df, gt_colname ,segment) {
  test_seg <- segment %in% df$vigor_Gene
  if (isFALSE(test_seg)) {
    final_gtype <- order$gt[order$vigor_Gene == segment]
  } else {
    idx <- which(df$vigor_Gene == segment)
    sub_seg <- df[idx,]
    sub_seg <- sub_seg[mixedorder(sub_seg[[gt_colname]]),]
    
    final_gtype <- paste0(sub_seg[[gt_colname]], collapse = '')
  }
  return(final_gtype)
}

get_row_col_idx <- function(df, col_name, pattern, negate) {
  out_df <- data.frame()
  for (i in 1:length(col_name)) {
    x <- which(colnames(df) == col_name[i])
    y <- str_which(string = df[[x]], 
                   pattern = pattern, 
                   negate = negate)
    if ((length(x) > 0) & (length(y) > 0)) {
      idx_col <- vector()
      idx_row <- vector()
      for (j in 1:length(y)) {
        idx_col <- c(idx_col, x)
        idx_row <- c(idx_row, y[j])
      }
      sub_df <- data.frame(idx_col, idx_row)
      out_df <- rbind(out_df, sub_df)
    }
  }
  out_df$idx_row <- out_df$idx_row + 1
  return(out_df)
}

# Function to map column names to their positions and format strings
generate_col_strings <- function(in_df, in_col_names) {
  # Validate inputs
  if (!is.data.frame(in_df)){
    stop("The first argument must be a data frame.")
  }
  if (!is.character(in_col_names)) {
    stop("Column names must be provided as a character vector.")
  }
  
  # Find column positions
  col_positions <- match(in_col_names, colnames(in_df))
  
  # Handle missing columns
  if (any(is.na(col_positions))) {
    missing_cols <- in_col_names[is.na(col_positions)]
    warning("The following columns were not found in the dataframe: ",
            paste(missing_cols, collapse = ", "))
    # Remove missing columns from processing
    in_col_names <- in_col_names[!is.na(col_positions)]
    col_positions <- col_positions[!is.na(col_positions)]
  }
  
  # Generate formatted strings
  out_df <- in_df[in_col_names]
  for (i in 1:ncol(out_df)) {
    out_df[[i]] <- paste0("[", colnames(out_df)[i], "=", out_df[[i]], "]")
  }
  result <- apply(out_df, 1, paste, collapse = " ", simplify = T) %>% as.vector()
  return(result)
}

get_GC_content <- function(Sequence_in) {
  Number.of.Sequences <- length(Sequence_in)
  Sequence <- paste0(Sequence_in, collapse = '')
  split_Sequence <- str_split(string = Sequence, pattern = '', n = Inf) %>% unlist()
  freq <- count(split_Sequence)
  
  nG <- freq$freq[freq$x == 'G']
  nC <- freq$freq[freq$x == 'C']
  nA <- freq$freq[freq$x == 'A']
  nT <- freq$freq[freq$x == 'T']
  
  idx <- str_which(freq$x, '[^GCAT]')
  nN <- sum(freq$freq[idx])
  
  if (nN > 0) {
    sub_freq <- freq[idx,]
    N <- paste0(sub_freq$x, '(', sub_freq$freq, ')', collapse = ', ')
    #notes <- paste0('Sequence(s) contain ambiguous nucleotides: ', N)
  } else {
    N <- NA
    #notes <- 'No notes.'
  }
  
  nTotal <- nG + nC + nA + nT + nN
  GC <- (nG + nC) / nTotal * 100
  
  sub_df <- data.frame(Number.of.Sequences, nG,nC,nA,nT,nN,nTotal, GC, N)
  return(sub_df)
}

clear_and_create_dirs <- function(dir_list, create) {
  for (dir_path in dir_list) {
    if (dir.exists(dir_path)) {
      unlink(dir_path, 
             recursive = TRUE, force = TRUE)
    }
    
    if (isTRUE(create)) {
      dir.create(dir_path, 
                 recursive = TRUE, showWarnings = FALSE)
    }
    
  }
}

## OUT_DIR ----
out_dir <- paste0(in_dir, '/3_analysis_output')
clear_and_create_dirs(dir_list = out_dir, 
                      create = TRUE)

output_seq_dir <- paste0(out_dir, "/output_sequences")
dir.create(path = output_seq_dir, 
           showWarnings = FALSE)

table2asn_dir <- paste0(in_dir, "/4_table2asn_input")
dir.create(path = table2asn_dir, 
           showWarnings = FALSE)

## FILE_LIST ----
file <- list.files(path = in_dir, 
                   recursive = TRUE, full.names = TRUE)


## NOTE: Rerunning gb_submission_3.R may overwrite "SeqIDs" already submitted to GB 
## "Seq_IDs" are replaced by GB accessions
if (any(str_detect(string = file, pattern = '(A|a)ccessions'))) {
  file_name <- file[str_which(string = file, pattern = '(A|a)ccessions')]
  message("\t", "Accession file detected. Running this script may overwrite temporary 'SeqIDs' that will be linked to a final GenBank accession.")
  message("\t\t", "file_name: ", file_name)
  quit()
}

## JOB_FILE ----
in_job_file <- grep(pattern = '.+job_file.+\\.txt', 
                    x = file, value = TRUE)
in_job_file <- grep(pattern = 'bowtie', 
                    invert = TRUE, 
                    x = in_job_file, value = TRUE)

job_file <- read.table(file = in_job_file, 
                       header = FALSE, sep = '\t', 
                       col.names = c('in_fasta', 'Read_1', 'Read_2', 'OUT_DIR'))

## INPUT ----
in_fasta <- job_file$in_fasta

in_sample_info <- grep(pattern = '.+analysis_input/Sample_Info.+\\.xlsx$',
                       x = file, value = TRUE)

in_sbt <- grep(pattern = '*.sbt$', 
               x = file, value = TRUE)

in_blast_1 <- grep(pattern = '.+blast_1/blast_.+\\.txt', 
                   x = file, value = TRUE)

in_blast_2 <- grep(pattern = '.+blast_2/blast_.+\\.txt', 
                   x = file, value = TRUE)

in_gene <- grep(pattern = '.+vigor/vigor4_.+\\.rpt_gene', 
                x = file, value = TRUE)

in_bowtie <- grep(pattern = '.+bowtie.+SAMTOOLS_COVERAGE.txt', 
                  x = file, value = TRUE)

in_depth <- grep(pattern = '.+bowtie.+SAMTOOLS_MEAN_DEPTH.txt', 
                 x = file, value = TRUE) 

in_seqkit <- grep(pattern = '.+bowtie.+SEQKIT_STATS.txt', 
                  x = file, value = TRUE) 

## FASTA ----
# pattern <- paste0(project_name, ".+$")
# file_string <- str_extract(string = in_fasta, pattern = pattern)
# file_string <- paste0("\t\t", file_string , collapse = "\n")
message("\n\t", "Reading input assembly sequences (IN_FASTA)...")

fasta <- concat_fas(file_list = in_fasta, 
                    add_file_name = TRUE, full_name = TRUE)

colnames(fasta) <- c('seq.ID', 'Sequence', 'in_fasta')

o_headers <- read.xlsx(xlsxFile = in_sample_info, 
                       sheet = "Original_Headers")

fasta <- merge(x = fasta, y = o_headers, 
               by = "seq.ID", all.x = TRUE, sort = FALSE)

fasta$Assembly_ID <- file_path_sans_ext(fasta$in_fasta) %>% basename()

fasta$Length_SeqID_NT <- str_length(string = fasta$Sequence)

## VIGOR gene ----
# pattern <- paste0(project_name, ".+$")
# file_string <- str_extract(string = in_gene, pattern = pattern)
# file_string <- paste0("\t\t", file_string , collapse = "\n")
message("\n\t", "Reading input VIGOR results (IN_GENE)...")

gene <- concat_tables(file_list = in_gene, header = TRUE)

#RVX
gene$RVX <- str_remove_all(gene$Ref_DB, '^rtv|_db$') %>% str_to_upper()
gene$RVX[gene$RVX == ""] <- "X"
gene$RVX <- paste0('RV', gene$RVX)

#organism
gene$organism <- str_replace(string = gene$RVX, 
                             pattern = '^RV', 
                             replacement = 'Rotavirus ')

colnames(gene) <- paste0('vigor_', colnames(gene))
names(gene)[names(gene) == 'vigor_Seq_ID'] <- 'seq.ID'
names(gene)[names(gene) == 'vigor_Notes'] <- 'Notes'


## GENE_3 ----
gene$gene_3 <- gene$vigor_Gene
idx <- which(gene$gene_3 == 'NSP5'|gene$gene_3 == 'NSP6')
seq_id <- unique(gene$seq.ID[idx])

for (i in 1:length(seq_id)) {
  idx <- which(gene$seq.ID == seq_id[i])
  gene$gene_3[idx] <- paste0('NSP', paste0(str_remove(mixedsort(unique(gene$gene_3[idx])), 'NSP'), collapse = '-'))
}

## LEN_ORF, ORF.2, & P_COV ----
## link NSP6 to corresponding NSP5 values for genotyping
len_orf <- vector()
orf.2 <- vector()
p_cov <- vector()

for (i in 1:length(gene$vigor_Gene)){
  len_orf[i] <- gene$vigor_Length_ORF_NT[i]
  orf.2[i] <- gene$vigor_ORF[i]
  p_cov[i] <- gene$vigor_Percent_Coverage[i]
}

gene <- data.frame(gene, len_orf, orf.2, p_cov, stringsAsFactors = FALSE)

## BLAST_1 ----
# pattern <- paste0(project_name, ".+$")
# file_string <- str_extract(string = in_blast_1, pattern = pattern)
# file_string <- paste0("\t\t", file_string , collapse = "\n")
message("\n\t", "Reading input BLAST_1-NCBI nt results (IN_BLAST_1)...")

blast_1 <- concat_tables(in_blast_1, header = TRUE)
colnames(blast_1) <- paste0('blast_1_', colnames(blast_1))
names(blast_1)[names(blast_1) == 'blast_1_qseqid'] <- 'seq.ID'

## TOP_5 HITS (BLAST_1)
x <- 5
top_x <- get_top_x(df = blast_1, x = x, grp_by_col = 'seq.ID')

## TOP_1 HITS (BLAST_1)
blast_1 <- get_top_x(df = blast_1, x = 1, grp_by_col = 'seq.ID')

## SACC (BLAST_1 TOP_5 HITS)
sacc <- get_sacc(df = top_x, sacc_col = 'blast_1_sacc')
out_sacc <- paste0(out_dir, '/blast_1_top_', x, '_sacc.txt')

write.table(x = sacc, file = out_sacc, sep = '\t', row.names = FALSE, 
            col.names = FALSE, quote = FALSE)

## BLAST_2 ----
# pattern <- paste0(project_name, ".+$")
# file_string <- str_extract(string = in_blast_2, pattern = pattern)
# file_string <- paste0("\t\t", file_string , collapse = "\n")
message("\n\t", "Reading input BLAST_2-RVA results (IN_BLAST_2)...")

blast_2 <- concat_tables(in_blast_2, header = T)
colnames(blast_2) <- paste0('blast_2_', colnames(blast_2))
names(blast_2)[names(blast_2) == 'blast_2_qseqid'] <- 'seq.ID'

freq <- count(df = blast_2, 
              vars = c("seq.ID", "blast_2_genotype"))

freq <- freq[mixedorder(freq$seq.ID),]
rownames(freq) <- NULL

## blast_2 hits (genotype frequencies)
hits <- data.frame()
u_seq.ID <- unique(blast_2$seq.ID)
for (i in 1:length(u_seq.ID)) {
  idx <- which(freq$seq.ID == u_seq.ID[i])
  sub_freq <- freq[idx,]
  blast_2_freq <- paste0(sub_freq$blast_2_genotype, '(', sub_freq$freq, ')', collapse = ', ')
  if (length(sub_freq$blast_2_genotype) > 1) {
    blast_2_note <- paste0("Variable blast_2_genotype hits (", blast_2_freq, ')')
  } else {
    blast_2_note <- 'No notes.'
  }
  df <- data.frame('seq.ID' = u_seq.ID[i], blast_2_freq, blast_2_note)
  hits <- rbind(hits, df)
}

## blast_2 tophit only
blast_2 <- get_top_x(df = blast_2, x = 1, grp_by_col = 'seq.ID')
blast_2 <- merge_dfs(x = blast_2, y = hits, by = 'seq.ID', all = TRUE)

## BOWTIE ----
if (length(in_bowtie) > 0) {
  df <- fasta[c('Assembly_ID','seq.ID')]
  if (length(in_bowtie) > 0) {
    ## BOWTIE in_bowtie ----
    # pattern <- paste0(project_name, ".+$")
    # file_string <- str_extract(string = in_bowtie, pattern = pattern)
    # file_string <- paste0("\t\t", file_string , collapse = "\n")
    message("\n\t", "Reading input Bowtie/Samtools coverage results (IN_BOWTIE)...")

    bowtie <- concat_tables(in_bowtie, header = TRUE)
    colnames(bowtie) <- paste0('bowtie_', colnames(bowtie))
    names(bowtie)[names(bowtie) == 'bowtie_rname'] <- 'seq.ID'
    
    bowtie <- merge_dfs(x = df, y = bowtie, by = 'seq.ID', all = TRUE)
    
    ## BOWTIE in_depth ----
    # pattern <- paste0(project_name, ".+$")
    # file_string <- str_extract(string = in_depth, pattern = pattern)
    # file_string <- paste0("\t\t", file_string , collapse = "\n")
    message("\n\t", "Reading input Bowtie/Samtools depth results (IN_DEPTH)...")
    
    depth <- concat_tables(in_depth, header = TRUE)
    colnames(depth) <- paste0('bowtie_', colnames(depth))
    names(depth)[names(depth) == 'bowtie_BASE'] <- 'Assembly_ID'
    
    ## BOWTIE in_seqkit ----
    # pattern <- paste0(project_name, ".+$")
    # file_string <- str_extract(string = in_seqkit, pattern = pattern)
    # file_string <- paste0("\t\t", file_string , collapse = "\n")
    message("\n\t", "Reading input Bowtie/SeqKit stats (IN_SEQKIT)...")
    
    seqkit <- concat_tables(in_seqkit, header = TRUE)
    colnames(seqkit) <- paste0('seqkit_', colnames(seqkit))
  }
}


## RVA Cutoffs (RCWG) ----
cutoff <- data.frame("vigor_Gene" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", 
                                      "NSP1", "NSP2", "NSP3", "NSP4", "NSP5"), 
                     "final.genotype" = c("Gx", "P[x]", "Ix", "Rx", "Cx", "Mx", 
                                          "Ax", "Nx", "Tx", "Ex", "Hx"), 
                     "cutoff" = c(80, 80, 85, 83, 84, 81, 79, 85, 85, 85, 91), 
                     stringsAsFactors = FALSE)


########## TABLE #########
## merge: gene, blast_1, blast_2, cutoff by 'Header'; this 'table' will serve as base for raw/unedited input values (including NAs -> no vigor/blast produced)

table <- merge(x = cutoff, y = gene, 
               by = "vigor_Gene", all = TRUE, sort = FALSE)
table <- merge(x = table, y = blast_1, 
               by = 'seq.ID', all = TRUE, sort = FALSE)
table <- merge(x = table, y = blast_2, 
               by = 'seq.ID', all = TRUE, sort = FALSE)
table <- merge(x = fasta, y = table, 
               by = 'seq.ID', all = TRUE, sort = FALSE)
if (length(in_bowtie) > 0) {
  table <- merge(x = bowtie, y = table, 
                 by = c('seq.ID', 'Assembly_ID'), 
                 all = TRUE, sort = FALSE)
}



########## SUMMARY ########
## extract stats from 'table' which produced vigor results ("gene.ID")
idx <- which(is.na(table$vigor_ORF_ID) == FALSE)   
summary <- table[idx,]

idx <- which(is.na(summary$Notes))
summary$Notes[idx] <- ""

idx <- which(summary$bowtie_meandepth < 30)

summary$Notes[idx] <- paste0(summary$Notes[idx], 
                             "; Average Read Depth (", 
                             summary$bowtie_meandepth[idx], 
                             ") < 30x")


## GENOTYPING ----
## determine if 'blast_2' can be assigned as final 'genotype' and 
## whether it will be 'included' in final 'constellation'
message("\n\t", "Determining RVA genotypes...")
order <- data.frame("vigor_Gene" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", 
                                     "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6"),
                    "gt_na" = c("Gna", "P[na]", "Ina", "Rna", "Cna", "Mna", 
                                "Ana", "Nna", "Tna", "Ena", "Hna", "Hna"), 
                    "gt_nt" = c("Gnt", "P[nt]", "Int", "Rnt", "Cnt", "Mnt", 
                                "Ant", "Nnt", "Tnt", "Ent", "Hnt", "Hnt"), 
                    stringsAsFactors = FALSE)

included <- vector()
for (i in 1:length(summary$vigor_ORF_ID)) {
  ## Skip NSP6
  if (summary$vigor_Gene[i] == 'NSP6') {
    included[i] = NA
    next(i)
  }
  ## NO blast_2 hit
  if (is.na(summary$blast_2_genotype[i])) {
    included[i] = "no"
    summary$note[i] <- paste0(summary$note[i], "; No blast_2 hit")
    idx <- which(order$vigor_Gene == summary$vigor_Gene[i])
    summary$final.genotype[i] <- order$gt_na[idx]
    summary$blast_2_genotype[i] <- order$gt_na[idx]
    next(i)   ### continue to next iteration of i
  }
  
  ## Variable blast_2_genotype hits
  if (isTRUE(str_detect(summary$blast_2_note[i], 'Variable'))) {
    included[i] = "yes"
    summary$note[i] <- paste0(summary$note[i], '; ',summary$blast_2_note[i])
  }
  
  ## Non-RVA
  if (summary$vigor_RVX[i] != 'RVA') {
    included[i] = "no"
    summary$Notes[i] <- paste0(summary$Notes[i], '; ', 'Non-RVA (vigor_RVX =', summary$vigor_RVX[i], ')')
    idx <- which(order$vigor_Gene == summary$vigor_Gene[i])
    summary$final.genotype[i] <- order$gt_nt[idx]
    next(i)
  }
  
  ## Gap or frame shift
  if (isTRUE(str_detect(summary$vigor_Location[i], pattern = ","))) {
    included[i] = "no"
    idx <- which(order$vigor_Gene == summary$vigor_Gene[i])
    summary$final.genotype[i] <- order$gt_nt[idx]
    next(i)
  }
  
  ## GENOTYPING complete ORF----
  if (summary$orf.2[i] == 'complete') {
    ## if blast_2_pi >= cutoff, assign blast_2_genotype as genotype
    if (summary$blast_2_pident[i] >= summary$cutoff[i]) {   
      included[i] = "yes"
      summary$final.genotype[i] <- summary$blast_2_genotype[i]
    } else {
      ## if blast_2_pi < cutoff, assign 'NT' as genotype
      included[i] = "yes"
      idx <- which(order$vigor_Gene == summary$vigor_Gene[i])
      summary$final.genotype[i] <- order$gt_nt[idx]
      summary$Notes[i] <- paste0(summary$Notes[i], "; Potential novel genotype will require further review by RCWG for confirmation; ",
                                'RCWG criteria for RVA genotype classification of complete ORFs: (1) complete ORF including start/stop codon, and (2) percent identity >= RCWG cutoff. This sequence contains a complete ORF, but cannot be classified because ',
                                "percent identity (", summary$blast_2_pident[i], ") ",
                                "< RCWG ", summary$blast_2_segment[i], " cutoff ", 
                                "(", summary$cutoff[i], ")")
    }
    ## GENOTYPING partial ORF ----
  } else {
    test_len <- summary$len_orf[i] >= 500
    test_pcov <- summary$p_cov[i] >= 50
    
    old_cutoff <- summary$cutoff[i]
    summary$cutoff[i] <- summary$cutoff[i] + 2
    test_pi <- summary$blast_2_pident[i] >= summary$cutoff[i]
    
    ## if all are TRUE, assign blast_2_genotype
    if (test_len & test_pcov & test_pi) {
      included[i] = "yes"
      summary$final.genotype[i] <- summary$blast_2_genotype[i]
      
    } else {
      included[i] = "no"
      idx <- which(order$vigor_Gene == summary$vigor_Gene[i])
      summary$final.genotype[i] <- order$gt_nt[idx]
      
      n = 0
      note <- vector()
      if (test_len == F) {
        note <- c(note, paste0("ORF length (", summary$len_orf[i], " nt) < 500 nt"))
      }
      if (test_pcov == F) {
        note <- c(note, paste0("ORF coverage (", summary$p_cov[i], "%) < 50%"))
      }
      if (test_pi == F) {
        note <- c(note, paste0("percent identity (", summary$blast_2_pident[i], "%) < RCWG ", summary$blast_2_segment[i], " cutoff + 2% (", summary$cutoff[i], "%)"))
      }
      note <- paste0(note, collapse = ', ')
      #note <- paste0('RCWG criteria for RVA genotype classification of partial ORFs: (1) ORF length >= 500 nt, (2) ORF coverage >= 50%, and (3) percent identity >= RCWG Cutoff + 2%. This sequence contains a partial ORF and cannot be classified because ', note)
      note <- paste0(summary$vigor_Gene[i], " (", summary$vigor_ORF_ID[i], ') contains a partial ORF and cannot be classified because ', note)
      
      summary$Notes[i] <- paste0("; ", note)
    }
  }
}

summary$Notes <- str_remove(string = summary$Notes,
                            pattern = '^; ')

## genotype without brackets
summary$gtype <- gsub(pattern = "\\[|\\]", 
                      replacement = "", 
                      x = summary$final.genotype)

## add 'included' to summary
summary <- data.frame(summary, 
                      included, stringsAsFactors = FALSE)

## Do not include seqs with multiple ORFs that could not be genotyped
idx_1 <- str_which(summary$vigor_ORF_ID, '[a-z]$')
idx_2 <- str_which(summary$final.genotype, 'x')

idx <- intersect(idx_1, idx_2)

if (length(idx) > 0) {
  summary$included[idx] <- 'no'
}


## C_TABLE -----
message("\n\t", "Generating genotype constellations (C_TABLE)...")
c_table <- data.frame()
u_Assembly_ID <- mixedsort(unique(summary$Assembly_ID))

for (i in 1:length(u_Assembly_ID)) {
  Assembly_ID = u_Assembly_ID[i]
  
  idx <- which(summary$Assembly_ID == Assembly_ID)
  sub_summary <- summary[idx,]
 
  sub_summary <- split(x = sub_summary, 
                       f = sub_summary$vigor_RVX)
  sub_c_table <- data.frame()
  
  for (j in 1:length(sub_summary)) {
    order <- data.frame("vigor_Gene" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6"),
                        "gt" = c("Gna", "P[na]", "Ina", "Rna", "Cna", "Mna", "Ana", "Nna", "Tna", "Ena", "Hna", "Hna"), stringsAsFactors = F)
    sub <- sub_summary[[j]]
    
    ## constellation_NT ----
    VP7 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'VP7')
    VP4 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'VP4')
    VP6 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'VP6')
    VP1 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'VP1')
    VP2 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype',segment = 'VP2')
    VP3 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'VP3')
    NSP1 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'NSP1')
    NSP2 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'NSP2')
    NSP3 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'NSP3')
    NSP4 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'NSP4')
    NSP5 <- get_gtype_simple(df = sub, gt_colname = 'final.genotype', segment = 'NSP5')
    
    constellation_NT <- c(VP7, VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, NSP4, NSP5)
    constellation_NT <- paste0(constellation_NT, collapse = '-') 
    
    ## constellation_blast2 ----
    VP7 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'VP7')
    VP4 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'VP4')
    VP6 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'VP6')
    VP1 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'VP1')
    VP2 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype',segment = 'VP2')
    VP3 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'VP3')
    NSP1 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'NSP1')
    NSP2 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'NSP2')
    NSP3 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'NSP3')
    NSP4 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'NSP4')
    NSP5 <- get_gtype_simple(df = sub, gt_colname = 'blast_2_genotype', segment = 'NSP5')
    
    constellation_blast_2 <- c(VP7, VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, NSP4, NSP5)
    constellation_blast_2 <- paste0(constellation_blast_2, collapse = '-')
    
    order <- data.frame("vigor_Gene" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6"),
                        "gt" = c("Gx", "P[x]", "Ix", "Rx", "Cx", "Mx", "Ax", "Nx", "Tx", "Ex", "Hx", "Hx"), stringsAsFactors = F)
    sub <- sub_summary[[j]]
    
    idx_yes <- which((sub$included == 'yes') | (is.na(sub$included)))
    
    sub <- sub[idx_yes,]
    
    vigor_RVX <- names(sub_summary)[j]
    
    analysis_output <- length(unique(sub$seq.ID))
    
    ORF.Complete <- length(which(sub$vigor_ORF == 'complete'))
    ORF.Partial <- length(which(sub$vigor_ORF == 'partial'))
    
    Missing <- setdiff(order$vigor_Gene, sub$vigor_Gene) %>% 
      paste0(collapse = ', ') %>% na_if(y = "")
    
    ## constellation_final.genotype ----
    VP7 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'VP7')
    VP4 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'VP4')
    VP6 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'VP6')
    VP1 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'VP1')
    VP2 <- get_gtypes(df = sub, gt_colname = 'final.genotype',segment = 'VP2')
    VP3 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'VP3')
    NSP1 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'NSP1')
    NSP2 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'NSP2')
    NSP3 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'NSP3')
    NSP4 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'NSP4')
    NSP5 <- get_gtypes(df = sub, gt_colname = 'final.genotype', segment = 'NSP5')
    
    constellation_final.genotype <- c(VP7, VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, NSP4, NSP5)
    constellation_final.genotype <- str_remove_all(string = constellation_final.genotype, 
                                                   pattern = '\\**\\([[:alnum:],:\\*]+\\)|\\*') %>% unique()
    constellation_final.genotype <- paste0(constellation_final.genotype, collapse = '-')  
    
    GxPx <- paste0(unique(c(VP7, VP4)), collapse = '') %>% str_remove_all('\\**\\([[:alnum:],:\\*]+\\)|\\*')
    
    c_row <- data.frame(GxPx,
                        VP7, VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, NSP4, NSP5, 
                        constellation_final.genotype,
                        constellation_NT, 
                        constellation_blast_2, 
                        Assembly_ID, 
                        vigor_RVX,
                        analysis_output, 
                        ORF.Complete, 
                        ORF.Partial, 
                        Missing)
    
    sub_c_table <- rbind(sub_c_table, 
                         c_row)
  }
  c_table <- rbind(c_table, 
                   sub_c_table)
}

## for NSP6 assign same genotype as NSP5
idx <- which(summary$vigor_Gene == 'NSP6')
if (length(idx)) {
  for (i in 1:length(idx)) {
    seq.id <- summary$seq.ID[idx[i]]
    idx_nsp5 <- which((summary$seq.ID == seq.id) & (summary$vigor_Gene == 'NSP5'))
    summary$final.genotype[idx[i]] <- paste0(mixedsort(unique(summary$final.genotype[idx_nsp5])), collapse = '')
    
    if (any(str_detect(summary$included[idx_nsp5], 'yes'))) {
      summary$included[idx[i]] <- 'yes'
    } else {
      summary$included[idx[i]] <- 'no'
    }
  }
}


## SAMPLE_INFO ----
message("\n\t", "Reading sample info table (IN_SAMPLE_INFO)...")
sample_info <- read.xlsx(xlsxFile = in_sample_info, 
                         sheet = 'Assembly_Info', 
                         colNames = TRUE, na.strings = '', detectDates = FALSE)

## Merge c_table and sample_info ----
col_names <- union(names(c_table), names(sample_info))
c_table <- merge(x = c_table, 
                 y = sample_info, 
                 by = 'Assembly_ID', all = TRUE, sort = FALSE)[,col_names]

idx <- which(is.na(c_table$analysis_output))
c_table$analysis_output[idx] <- 0

## Check if any prefix were provided for common_name 
idx <- which(!is.na(c_table$common_name.prefix))
if (length(idx) > 0) {
  c_table$common_name <- paste0(c_table$common_name.prefix, '-', 
                                c_table$NGS_ISOLATE)
} else {
  c_table$common_name <- c_table$NGS_ISOLATE
}


## STRAIN ----
message("\n\t", "Generating rotavirus 'strain' values...")
idx <- which(is.na(c_table$vigor_RVX))
c_table$vigor_RVX[idx] <- "RVX"

idx <- which(is.na(c_table$GxPx))
c_table$GxPx[idx] <- "GxP[x]"

strain <- paste0(c_table$vigor_RVX, '/', 
                 c_table$host.common, '-', c_table$type, '/',
                 c_table$iso_alpha3_code, '/', 
                 c_table$common_name, '/', 
                 str_extract(c_table$collection.date, '[0-9X]{4}'), '/', 
                 c_table$GxPx)

c_table <- data.frame(strain, c_table, stringsAsFactors = F)



## Merge summary and c_table by 'ngs_sample_ID' ----
summary <- merge(x = summary, y = c_table, 
                 by = c('Assembly_ID', 'vigor_RVX'), 
                 sort = FALSE)

## FINAL HEADERS ----
message("\n\t", "Generating 'final.headers' (GenBank)...")
## FINAL HEADERS p_table ----
p_table <- data.frame("gene_3" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", 
                                   "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", 
                                   "NSP5-6", "NSP6"), 
                      "product_2" = c("capsid glycoprotein VP7 (VP7) gene", 
                                      "outer capsid spike protein VP4 (VP4) gene", 
                                      "inner capsid protein VP6 (VP6) gene", 
                                      "RNA-dependent RNA polymerase VP1 (VP1) gene", 
                                      "core capsid protein VP2 (VP2) gene", 
                                      "RNA capping protein VP3 (VP3) gene", 
                                      "non-structural protein 1 (NSP1) gene",
                                      "non-structural protein 2 (NSP2) gene",
                                      "non-structural protein 3 (NSP3) gene", 
                                      "non-structural protein 4 (NSP4) gene", 
                                      "non-structural protein 5 (NSP5) gene",
                                      "non-structural protein 5 (NSP5) and non-structural protein 6 (NSP6) genes",
                                      "non-structural protein 6 (NSP6) gene"),
                      "segment" = c('9','4','6', '1','2','3','5','8','7','10','11', '11', '11'), stringsAsFactors = F)

## merge summary and p_table by 'gene_2' (adds 'product_2' and 'segment' to summary)
col_names <- union(names(summary), names(p_table))
summary <- merge(x = summary, y = p_table, 
                 by = 'gene_3', all.x = TRUE, sort = FALSE)[,col_names]


## FINAL HEADERS Definition ----
summary$final.header <- paste0(summary$vigor_organism, 
                               ' strain ', summary$strain,
                               ' ', summary$product_2,
                               ', ' ,summary$orf.2, ' cds')

## FINAL HEADERS Definition with NSP6 ----
idx <- which(summary$vigor_Gene == 'NSP6')
if (length(idx)) {
  for (i in 1:length(idx)) {
    seq.id <- summary$seq.ID[idx[i]]
    idx_nsp5 <- which((summary$seq.ID == seq.id) & (summary$vigor_Gene == 'NSP5'))
    
    if (any(str_detect(summary$included[idx_nsp5], 'yes'))) {
      idx_yes <- which(summary$included[idx_nsp5] == 'yes')
      summary$final.header[idx[i]] <- summary$final.header[idx_nsp5[idx_yes]]
    } else {
      summary$final.header[idx[i]] <- paste0(summary$vigor_organism[idx[i]], 
                                             ' strain ', summary$strain[idx[i]],
                                             ' ', "non-structural protein 6 (NSP6) gene", 
                                             ', ' ,summary$orf.2[idx[i]], ' cds')
    }
  }
}

## SUMMARY_FULL ----
## same as 'summary' except its sorted by isolate, segment, final.genotype
summary_full <- data.frame()
order <- data.frame("gene_3" = c("VP7", "VP4", "VP6", "VP1", "VP2", "VP3", 
                                 "NSP1", "NSP2", "NSP3", "NSP4", "NSP5","NSP5-6"), 
                    stringsAsFactors = FALSE)

u_Assembly_ID <- unique(c_table$Assembly_ID)
for (i in 1:length(u_Assembly_ID)) {
  Assembly_ID <- u_Assembly_ID[i]
  idx <- which(summary$Assembly_ID == Assembly_ID)
  sub_strain <- summary[idx,]
  for (j in 1:length(order$gene_3)){
    idx <- which(sub_strain$gene_3 == order$gene_3[j])
    sub_seg <- sub_strain[idx,]
    sub_seg <- sub_seg[mixedorder(sub_seg$final.genotype),]
    sub_seg <- sub_seg[mixedorder(sub_seg$vigor_ORF_ID),]
    summary_full <- rbind(summary_full, sub_seg, stringsAsFactors = FALSE)
  }
}

## SEQID (final) ----
message("\n\t", "Generating final 'SeqID's...")
seqid <- data.frame('seq.ID' = unique(summary_full$seq.ID), 
                    stringsAsFactors = FALSE)
max_len <- rownames(seqid) %>% str_length() %>% max()
seqid$SeqID <- str_pad(string = rownames(seqid), 
                       width = max_len, 
                       side = 'left', pad = '0')
seqid$SeqID <- paste0('Seq', seqid$SeqID)

SeqID <- vector()
for (i in 1:length(summary_full$seq.ID)) {
  idx <- which(seqid$seq.ID == summary_full$seq.ID[i])
  SeqID[i] <- seqid$SeqID[idx]
}
summary_full<- data.frame(SeqID, summary_full, 
                          stringsAsFactors = FALSE)

## add SeqID to definition
summary_full$final.header.orf <- paste0(summary_full$SeqID,'_', 
                                        summary_full$vigor_ORF_ID ,' ', 
                                        summary_full$final.header)
summary_full$final.header <- paste0(summary_full$SeqID, ' ', 
                                    summary_full$final.header)

## GC_SEQUENCE ----
message("\n\t", "Calculating %GC by sequence (GC_SEQUENCE)...")
GC_SeqID <- data.frame('Assembly_ID' = summary_full$Assembly_ID, 
                       'SeqID' = summary_full$SeqID,
                       'seq.ID' = summary_full$seq.ID,
                       'strain' = summary_full$strain,
                       lapply(summary_full$Sequence, FUN = get_GC_content) %>% 
                         bind_rows())
summary_full$GC_SeqID <- GC_SeqID$GC

GC_SeqID <- GC_SeqID[which(!duplicated(GC_SeqID$SeqID)),]

## GC_ASSEMBLY ----
message("\n\t", "Calculating %GC by assembly (GC_ASSEMBLY)...")
split_summary_full <- summary_full[which(!duplicated(summary_full$SeqID)),]
split_summary_full <- split(x = split_summary_full, 
                            f = split_summary_full$Assembly_ID)

GC_Assembly <- data.frame('Assembly_ID' = lapply(split_summary_full, "[[", "Assembly_ID") %>%
                            lapply(FUN = unique) %>% 
                            unlist(use.names = FALSE),
                          sapply(split_summary_full,
                                function(x) apply(x[c("Sequence")], 2, get_GC_content)) %>%
                            bind_rows())

c_table <- merge(x = c_table, y = GC_Assembly[c('Assembly_ID', 'GC', 'nTotal')], 
                 by = 'Assembly_ID', all = TRUE, sort = FALSE)

names(c_table)[names(c_table) == 'GC'] <- 'GC_Assembly'
names(c_table)[names(c_table) == 'nTotal'] <- 'Total.Length.NT'

## SUMMARY_OUT ----
summary_out <- summary_full

## TABLE2ASN ----
message("\n\t", "Generating TABLE2ASN input files...")
base_name <- paste0(table2asn_dir, '/' , project_name)


## ## TABLE2ASN df ----
idx_1 <- which(summary_full$vigor_Gene != 'NSP6')

idx_2 <- which(summary_full$included != 'no')

#idx <- intersect(x = idx_1, idx_2)
idx <- unique(idx_1, idx_2)
df <- summary_full[idx,]

## ## TABLE2ASN .fsa ----
Header <- df$final.header
Sequence <- df$Sequence
fas <- build_fas(Header = Header, Sequence = Sequence)
out_fas <- paste0(base_name, '.fsa')
writeFasta(fdta = fas, out.file = out_fas)

## TABLE2ASN .fsa ORF_NT ----

Header <- df$final.header.orf
Sequence <- df$vigor_ORF_NT
fas <- build_fas(Header = Header, Sequence = Sequence)
out_fas <- paste0(base_name, '_ORF_NT.fas')
writeFasta(fdta = fas, out.file = out_fas)

## TABLE2ASN .fsa ORF_AA ----
Header <- df$final.header.orf
Sequence <- df$vigor_ORF_AA
fas <- build_fas(Header = Header, Sequence = Sequence)
out_fas <- paste0(base_name, '_ORF_AA.fas')
writeFasta(fdta = fas, out.file = out_fas)

## OUTPUT SEQUENCES By Segment ----
seq_dir <- paste0(output_seq_dir, "/by_segment")
dir.create(path = seq_dir, 
           showWarnings = FALSE, recursive = TRUE)

df$run_gene <- paste0(df$NGS_RUN_ID, "__", df$NGS_METHOD_ID, "__", df$vigor_Gene)
by_name <- "run_gene"
u_name <- df[[by_name]] %>% unique()

for(i in 1:length(u_name)) {
  k <- which(df[[by_name]] == u_name[i])
  sub_df <- df[k,]
  
  ## Full-length
  in_col_names <- c("NGS_RUN_ID", "NGS_METHOD_ID", "Assembly_ID", "SeqID", 
                    "seq.ID", "vigor_ORF_ID", "vigor_Gene", "vigor_ORF", 
                    "final.genotype", "vigor_Location", "vigor_Codon_Start",
                    "constellation_final.genotype")
  
  header_info <- generate_col_strings(in_df = sub_df, 
                                      in_col_names = in_col_names)
  
  Header <- paste0(sub_df$Assembly_ID, "-", 
                   sub_df$SeqID, "-", 
                   sub_df$seq.ID, " ",
                   str_remove_all(sub_df$strain, pattern = "\\[|\\]"))
  
  Header <- paste0(Header, " ", header_info)
  Sequence <- sub_df$Sequence
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/', 
                    u_name[i], '__FULL-LENGTH.fas')
  writeFasta(fdta = fas, out.file = out_fas)
  
  ## ORF_NT
  header_info <- generate_col_strings(in_df = sub_df, 
                                      in_col_names = in_col_names)
  
  Header <- paste0(sub_df$Assembly_ID, "-", 
                   sub_df$SeqID, "-", 
                   sub_df$vigor_ORF_ID, " ",
                   str_remove_all(sub_df$strain, pattern = "\\[|\\]"))
  
  Header <- paste0(Header, " ", header_info)
  Sequence <- sub_df$vigor_ORF_NT
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/',
                    u_name[i], '__ORF_NT.fas')
  writeFasta(fdta = fas, out.file = out_fas)
  
  ## ORF_AA
  Sequence <- sub_df$vigor_ORF_AA
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/',
                    u_name[i], '__ORF_AA.fas')
  writeFasta(fdta = fas, out.file = out_fas)
}

## OUTPUT SEQUENCES By Sample ----
seq_dir <- paste0(output_seq_dir, "/by_sample")
dir.create(path = seq_dir, 
           showWarnings = FALSE, recursive = TRUE)

df$run_sample <- paste0(df$NGS_RUN_ID, "__", 
                        df$NGS_METHOD_ID, "__", 
                        df$Assembly_ID)
by_name <- "run_sample"
u_name <- df[[by_name]] %>% unique()

for(i in 1:length(u_name)) {
  k <- which(df[[by_name]] == u_name[i])
  sub_df <- df[k,]
  
  ## Full-length
  in_col_names <- c("NGS_RUN_ID", "NGS_METHOD_ID", "Assembly_ID", "SeqID", 
                    "seq.ID", "vigor_ORF_ID", "vigor_Gene", "vigor_ORF", 
                    "final.genotype", "vigor_Location", "vigor_Codon_Start",
                    "constellation_final.genotype")
  
  header_info <- generate_col_strings(in_df = sub_df, 
                                      in_col_names = in_col_names)
  
  Header <- paste0(sub_df$Assembly_ID, "-", 
                   sub_df$SeqID, "-", 
                   sub_df$seq.ID, " ",
                   str_remove_all(sub_df$strain, pattern = "\\[|\\]"))
  
  Header <- paste0(Header, " ", header_info)
  Sequence <- sub_df$Sequence
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/', 
                    u_name[i], '__FULL-LENGTH.fas')
  writeFasta(fdta = fas, out.file = out_fas)
  
  ## ORF_NT
  header_info <- generate_col_strings(in_df = sub_df, 
                                      in_col_names = in_col_names)
  
  Header <- paste0(sub_df$Assembly_ID, "-", 
                   sub_df$SeqID, "-", 
                   sub_df$vigor_ORF_ID, " ",
                   str_remove_all(sub_df$strain, pattern = "\\[|\\]"))
  
  Header <- paste0(Header, " ", header_info)
  Sequence <- sub_df$vigor_ORF_NT
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/',
                    u_name[i], '__ORF_NT.fas')
  writeFasta(fdta = fas, out.file = out_fas)
  
  ## ORF_AA
  Sequence <- sub_df$vigor_ORF_AA
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  out_fas <- paste0(seq_dir, '/',
                    u_name[i], '__ORF_AA.fas')
  writeFasta(fdta = fas, out.file = out_fas)
}

## TABLE2ASN .src ----
source <- df[c('SeqID', 'vigor_organism', 'NGS_ISOLATE', 'collection.date', 
                         'country', 'host','vigor_Gene', 'final.genotype','isolation.source',
                         'strain','constellation_final.genotype')]

source$constellation_final.genotype <- paste0('constellation: ', 
                                              source$constellation_final.genotype)

colnames(source) <- c('SeqID', 'organism', 'isolate', 'collection-date', 
                      'country', 'host', 'segment', 'genotype','isolation-source', 
                      'strain', 'note')

out_source <- paste0(base_name, '.src')
write.table(x = source, file = out_source, 
            append = FALSE, quote = FALSE, sep = '\t', col.names = TRUE, 
            row.names = FALSE)

## TABLE2ASN .tbl ----
summary_full$final.product <- str_remove(string = summary_full$vigor_Product, 
                                         pattern = 'putative ')
summary_full$final.product <- str_remove(string = summary_full$final.product, 
                                         pattern = ', .+$') ## N/C-terminal/fragment
summary_full$start <- str_extract(string = summary_full$vigor_Location, 
                                  pattern = '^<*[0-9]+')
summary_full$stop <- str_extract(string = summary_full$vigor_Location, 
                                 pattern = '>*[0-9]+$')

out_feature <- paste0(base_name, '.tbl')
if(file.exists(out_feature)) {
  unlink(out_feature)
}

SeqID <- unique(summary_full$SeqID)
idx <- which(is.na(summary_full$Notes))
summary_full$Notes[idx] <- ""

for (i in 1:length(SeqID)) {
  idx <- which(summary_full$SeqID == SeqID[i])
  sub <- summary_full[idx,]
  sink(file = out_feature,append = T)
  cat(paste0('>Features', ' ', SeqID[i], '\n'))
  for (j in 1:length(sub$vigor_Gene)) {
    if (sub$included[j] == 'no') {
      next(j)
    }
    cat(paste0(sub$start[j], '\t', sub$stop[j], '\t', 'gene', '\n'))
    cat(paste0('\t\t\t', 'gene', '\t', sub$vigor_Gene[j], '\n'))
  
    cat(paste0(sub$start[j], '\t', sub$stop[j], '\t', 'CDS', '\n'))
    cat(paste0('\t\t\t', 'codon_start', '\t', sub$vigor_Codon_Start[j], '\n'))
    #cat(paste0('\t\t\t', 'protein_id', '\t', sub$gene.ID[j], '\n'))
    cat(paste0('\t\t\t', 'gene', '\t', sub$vigor_Gene[j], '\n'))
    cat(paste0('\t\t\t', 'product', '\t', sub$final.product[j], '\n'))
  
    if (str_detect(sub$Notes[j], 'further review|No notes', negate = T)) {
      note = sub$Notes[j]
      cat(paste0('\t\t\t', 'note', '\t', note, '\n'))
    }
  }
  sink()
}

## SUMMARY OUT ----
## SUMMARY OUT Concatenate Notes ----
u_Assembly_ID <- c_table$Assembly_ID
Notes <- vector()

for (i in 1:length(u_Assembly_ID)){
  idx <- which(summary_out$Assembly_ID == u_Assembly_ID[i])
  df <- summary_out[idx,]
  idx <- which(df$Notes != '')
  if (length(idx) > 0) {
    out_note <- paste0(df$vigor_Gene[idx], 
                       "(", df$SeqID[idx], "/", df$vigor_ORF_ID[idx], "): ", 
                       df$Notes[idx])
    out_note <- paste0(out_note, collapse = ";; ")
  } else {
    out_note <- ""
  }
  Notes[i] <- out_note
}

c_table$Notes <- Notes

c_colnames <- c('NGS_INSTRUMENT_ID', 'NGS_RUN_ID', 'NGS_METHOD_ID', 
                'Assembly_N', 'Assembly_ID', 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 
                'NGS_SN', 'NGS_N', 'vigor_RVX', 'Notes',
                'original_input', 'analysis_input','analysis_output',
                'strain', 'GxPx',
                'VP7','VP4', 'VP6','VP1','VP2','VP3', 
                'NSP1','NSP2','NSP3','NSP4','NSP5',
                'constellation_final.genotype', 
                'constellation_NT', 
                'constellation_blast_2', 
                'ORF.Complete', 'ORF.Partial', 
                'Missing', 'GC_Assembly', 'Total.Length.NT')

## Keep all columns. Append behing c_colnames
# c_colnames <- c(c_colnames, 
#                 setdiff(x = colnames(c_table),
#                         y = c_colnames))

## SUMMARY OUT c_table ----
c_table <- c_table[c_colnames]

## SUMMARY OUT c_table get row and col idx  ----

col_name <- c('VP7','VP4','VP6','VP1','VP2','VP3',
             'NSP1','NSP2','NSP3','NSP4','NSP5')

pattern <- '\\*'
COMPLETE <- get_row_col_idx(df = c_table, col_name = col_name, 
                            pattern = pattern, negate = FALSE)

pattern <- '\\(0\\)|\\[*x\\]*\\([1-9]+\\)'
PARTIAL <- get_row_col_idx(df = c_table, col_name = col_name, 
                              pattern = pattern, negate = TRUE)

pattern <-  '\\(0\\)|\\[*x\\]*\\([1-9]+\\)'
MISSING <- get_row_col_idx(df = c_table, col_name = col_name, 
                          pattern = pattern, negate = FALSE)

for (i in 1:length(col_name)) {
  string <- c_table[[col_name[i]]]
  pattern <- '\\*|\\(0\\)'
  c_table[[col_name[i]]] <- str_remove_all(string = string, 
                                           pattern = pattern)
}


## WB ----
wb <- createWorkbook()
sheet_name <- 'CONSTELLATIONS'
x <- c_table

idx_col <- which(colnames(x) == 'Notes')

addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')

setColWidths(wb = wb, sheet = sheet_name, 
             cols = 2:ncol(x), widths = "auto")
setColWidths(wb = wb, sheet = sheet_name, 
             cols = idx_col, widths = 50)

ORANGE <- createStyle(fontColour = '#9C5700', fgFill = '#FFEB9C')
addStyle(wb = wb, sheet = sheet_name, 
         style = ORANGE, 
         rows = PARTIAL$idx_row, 
         cols = PARTIAL$idx_col, 
         gridExpand = FALSE, stack = TRUE)

RED <- createStyle(fontColour = '#9C0006', fgFill = '#FFC7CE')
addStyle(wb = wb, sheet = sheet_name, 
         style = RED, 
         rows = MISSING$idx_row, 
         cols = MISSING$idx_col, 
         gridExpand = FALSE, stack = TRUE)

GREEN <- createStyle(fontColour = '#006100', fgFill = '#C6EFCE')
addStyle(wb = wb, sheet = sheet_name, 
         style = GREEN, 
         rows = COMPLETE$idx_row, 
         cols = COMPLETE$idx_col, 
         gridExpand = FALSE, stack = TRUE)


## WB vigor ----
col_names <- setdiff(grep(pattern = 'vigor_', x = colnames(summary_out), value = T), 'vigor_ORF_ID')
col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_ID', 'Assembly_N', 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N','SeqID', 'seq.ID', 'vigor_ORF_ID', col_names)
vigor <- summary_out[col_names]

## WB blast_1 ----
col_names <- setdiff(grep(pattern = 'blast_1_', x = colnames(summary_out), value = T), 
                     c('blast_1_stitle', 'blast_1_sacc', 'blast_1_pident'))
col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_ID', 'Assembly_N', 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N','SeqID', 'seq.ID','blast_1_stitle', 'blast_1_sacc', 'blast_1_pident', col_names)
blast_1 <- summary_out[col_names]
blast_1 <- blast_1[!duplicated(blast_1$SeqID),]

## WB blast_2 ----
col_names <- setdiff(grep(pattern = 'blast_2_', x = colnames(summary_out), value = T), 
                     c('blast_2_def', 'blast_2_sacc', 'blast_2_pident', 'blast_2_segment', 'blast_2_genotype'))
col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_ID', 'Assembly_N', 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N','SeqID', 'seq.ID', 'blast_2_segment', 'blast_2_sacc', 'blast_2_genotype','blast_2_def', 'blast_2_pident', col_names)
blast_2 <- summary_out[col_names]
blast_2 <- blast_2[!duplicated(blast_2$SeqID),]

## WB bowtie ----
if (length(in_bowtie) > 0) {
  col_names <- grep(pattern = 'bowtie_', x = colnames(summary_out), value = T)
  col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_ID', 'Assembly_N', 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N','SeqID', 'seq.ID', col_names)
  bowtie <- summary_out[col_names]
  bowtie <- bowtie[!duplicated(bowtie$SeqID),]
}

if (length(in_bowtie) > 0) {
  
  col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_N', 'Assembly_ID', 
                 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N',
                 'SeqID', 'seq.ID','vigor_ORF_ID', 'Notes', 'included', 
                 'vigor_Gene','final.genotype', 'vigor_ORF', 'cutoff', 
                 'blast_2_genotype', 'blast_2_pident', 'Length_SeqID_NT', 
                 'vigor_Length_ORF_NT', 'vigor_Length_ORF_AA', 
                 'vigor_Percent_Coverage', 'vigor_Location', 'vigor_Product', 
                 'bowtie_numreads', 'bowtie_meandepth', 
                 'blast_1_sacc', 'blast_1_stitle', 
                 'blast_2_sacc', 'blast_2_def', 'final.header', 'original_header')
  summary_out <- summary_out[col_names]
} else {
  col_names <- c('NGS_RUN_ID', 'NGS_METHOD_ID', 'Assembly_N', 'Assembly_ID', 
                 'NGS_ISOLATE_SN', 'NGS_ISOLATE', 'NGS_SN', 'NGS_N',
                 'SeqID', 'seq.ID', 'vigor_ORF_ID', 'Notes', 'included', 
                 'vigor_Gene','final.genotype', 'vigor_ORF', 'cutoff', 
                 'blast_2_genotype', 'blast_2_pident', 'Length_SeqID_NT', 
                 'vigor_Length_ORF_NT', 'vigor_Length_ORF_AA', 
                 'vigor_Percent_Coverage', 'vigor_Location', 'vigor_Product',
                 'blast_1_sacc', 'blast_1_stitle', 
                 'blast_2_sacc', 'blast_2_def', 'final.header', 'original_header')
  
  summary_out <- summary_out[col_names]
}

## WB sheets ----
sheet_name <- 'SUMMARY'
x <- summary_out
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

idx_row <- which(x$Notes != '') + 1
idx_col <- which(colnames(x) == 'Notes')
align_text <- createStyle(halign = "left", valign = "top")
wrap_text <- createStyle(wrapText = TRUE)
yellow_bg <- createStyle(fgFill = '#FFEB9C')
red_fnt <- createStyle(fontColour = 'red')

addStyle(wb = wb, sheet = sheet_name, style = align_text, 
         rows = 1:nrow(x), cols = 1:ncol(x), gridExpand = TRUE, stack = TRUE)
addStyle(wb = wb, sheet = sheet_name, style = yellow_bg, 
         rows = idx_row, cols = 1:ncol(x), gridExpand = TRUE, stack = TRUE)
addStyle(wb = wb, sheet = sheet_name, style = red_fnt, 
         rows = idx_row, cols = idx_col, gridExpand = TRUE, stack = TRUE)
addStyle(wb = wb, sheet = sheet_name, style = wrap_text, 
         rows = 1:nrow(x), cols = idx_col, gridExpand = TRUE, stack = TRUE)
setColWidths(wb = wb, sheet = sheet_name, cols = idx_col, widths = 50)

sheet_name <- 'Sample_Info'
x <- sample_info
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

sheet_name <- 'VIGOR'
x <- vigor
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

sheet_name <- 'BLAST_1-NCBI_NT'
x <- blast_1
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

sheet_name <- 'BLAST_2-RVA'
x <- blast_2
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

if (length(in_bowtie) > 0) {
  sheet_name <- 'BOWTIE_Assembly'
  x <- depth
  addWorksheet(wb = wb, sheetName = sheet_name)
  writeDataTable(wb = wb, sheet = sheet_name, x = x, 
                 colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
  setColWidths(wb = wb, sheet = sheet_name, 
               cols = 1:ncol(x), widths = "auto")
  
  sheet_name <- 'BOWTIE_Sequence'
  x <- bowtie
  addWorksheet(wb = wb, sheetName = sheet_name)
  writeDataTable(wb = wb, sheet = sheet_name, x = x, 
                 colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
  setColWidths(wb = wb, sheet = sheet_name, 
               cols = 1:ncol(x), widths = "auto")
  
  sheet_name <- 'SEQKIT'
  x <- seqkit
  addWorksheet(wb = wb, sheetName = sheet_name)
  writeDataTable(wb = wb, sheet = sheet_name, x = x, 
                 colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
  setColWidths(wb = wb, sheet = sheet_name, 
               cols = 1:ncol(x), widths = "auto")
}

sheet_name <- 'GC_Assembly'
x <- GC_Assembly
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

sheet_name <- 'GC_Sequence'
x <- GC_SeqID
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")


out_summary <- paste0(in_dir,'/Summary_', basename(in_dir), '.xlsx')

message("\n\t", "Saving analysis summary report...", "\n",
        "\t\t", "OUT_SUMMARY: ", out_summary)

saveWorkbook(wb = wb, file = out_summary, overwrite = TRUE)
