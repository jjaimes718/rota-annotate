#!/usr/bin/Rscript

## PACKAGES  ----
required_pkgs <- c("stringr", "microseq", "openxlsx", "tools", "gtools", "optparse", "plyr")

missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
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

in_dir <- normalizePath(path = in_dir, mustWork = TRUE)

project_name <- basename(in_dir)

message('RUNNING: 2_analysis_input.R', 
        '\n\n\t', 'IN_DIR: ', in_dir)

## FUNCTIONS ----
concat_fas <- function(file_list, add_file_name, full_name){
  fas <- data.frame()
  for (i in 1:length(file_list)) {
    if (file.size(file_list[i]) == 0) {
      message("\n\t","concat_fas: following file is empty\n",
          paste0("\t\tfile_list[", i, "]: "), file_list[i], "\n")
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

rename_col <- function(df, old_name, new_name) {
  if (isTRUE(old_name %in% colnames(df))) {
    names(df)[names(df) == old_name] <- new_name
    return(df)
  } else {
    cat("rename_col: old_name='", old_name, "', not found in colnames.", sep = "")
    return(df)
  }
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

grep_na <- function(pattern, x, ...) {
  res <- grep(pattern, x, ...)
  if(length(res) == 0) NA else res
}

## OUTPUT DIR ----

analysis_in_dir <- paste0(in_dir, "/2_analysis_input")

message("\n\t", "ANALYSIS_IN_DIR: ", analysis_in_dir)

## FILE_LIST ----
file_list <- list.files(path = in_dir, 
                        full.names = TRUE, recursive = TRUE)

## O_INPUT ----
pattern <- "/2_analysis_input/original_input.+\\.txt$"

in_o_input <- grep(pattern = pattern, 
                   x = file_list, value = TRUE)

o_input <- read.table(file = in_o_input, 
                      header = TRUE, sep = "\t")

## SAMPLE_INFO ----
pattern <- "/2_analysis_input/Sample_Info.+\\.xlsx$"

in_info <- grep(pattern = pattern, 
                x = file_list, value = TRUE)

sheet_name <- "Sample_Info"

info <- read.xlsx(xlsxFile = in_info, 
                  sheet = sheet_name, 
                  colNames = TRUE, na.strings = "")

col_nm_info <- colnames(info)

## MERGE o_input and info ----
col_names <- union(names(o_input), 
                   names(info))

o_input <- merge(x = o_input, 
                 y = info, 
                 all = TRUE, sort = FALSE)[,col_names]

## O_INPUT collection.date ----
idx <- which(is.na(o_input$collection.date))

if (length(idx) > 0) {
  message("\n\t", "Missing Sample_Info 'collection.dates' values will be deduced from 'NGS_RUN_ID'")
  o_input$collection.date[idx] <- str_extract(string = o_input$NGS_RUN_ID[idx], 
                                              pattern = "^[0-9]{2}")
  o_input$collection.date[idx] <- paste0('20', o_input$collection.date[idx])
}


## O_INPUT country ----
idx <- which(is.na(o_input$country))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'country' values will be set to default 'United States of America'")
  o_input$country[idx] <- "United States of America"
}

## O_INPUT Iso_alpha3_code ----
idx <- which(is.na(o_input$iso_alpha3_code))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'iso_alpha3_code' values will be set to default 'USA'")
  o_input$iso_alpha3_code[idx] <- "USA"
}

## O_INPUT Host ----
idx <- which(is.na(o_input$host))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'host' values will be set to default 'Homo sapiens'")
  o_input$host[idx] <- "Homo sapiens"
}

## O_INPUT Host.common ----
idx <- which(is.na(o_input$host.common))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'host.common' values will be set to default 'Human'")
  o_input$host.common[idx] <- "Human"
}

## O_INPUT Type ----
idx <- which(is.na(o_input$type))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'type' values will be set to default 'wt'")
  o_input$type[idx] <- "wt"
}

## O_INPUT isolation.source ----
idx <- which(is.na(o_input$isolation.source))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'isolation.source' values will be set to default 'stool'")
  o_input$isolation.source[idx] <- "stool"
}

## O_INPUT Common_name.prefix ----
idx <- which(is.na(o_input$common_name.prefix))

if (length(idx) > 0) {
  message("\t", "Missing Sample_Info 'common_name.prefix' values will be set to default '' (empty string)")
  o_input$common_name.prefix[idx] <- ""
}


## FASTA in_fasta ----
in_fasta <- o_input$in_assembly %>% na.omit()

## FASTA fasta ----
fasta <- concat_fas(file_list = in_fasta, 
                    add_file_name = TRUE, full_name = TRUE)

## FASTA CTG_IDs ----
message("\n\t", "Generating FASTA Contig IDs...")

max_len <- rownames(fasta) %>% str_length() %>% max()

fasta$seq.ID <- str_pad(string = rownames(fasta), 
                        width = max_len, 
                        side = 'left', 
                        pad = '0')

fasta$seq.ID <- paste0('CTG_', fasta$seq.ID)

## FASTA original_header ----
fasta <- rename_col(df = fasta, 
                    old_name = "Header", 
                    new_name = "original_header")

## FASTA in_assembly ----
fasta$in_assembly <- fasta$file_name

## INPUT_ORIGINAL_SEQUENCES ----
freq <- count(df = fasta, vars = "in_assembly")

freq <- rename_col(df = freq, 
                   old_name = "freq", 
                   new_name = "original_input")

col_names <- union(names(o_input), 
                   names(freq))

o_input <- merge(x = o_input, 
                 y = freq, 
                 all = TRUE, sort = FALSE)[,col_names]

idx_na <- which(is.na(o_input$original_input))

o_input$original_input[idx_na] <- 0

## FASTA edit headers ----
message("\t", "Cleaning FASTA sequence headers...")
fasta$original_header <- str_remove_all(string = fasta$original_header, 
                                        pattern = "[\\(\\),\\[\\]]+")
## FASTA edit file_name ----
message("\t", "Cleaning FASTA file names...")
fasta$file_name <- str_remove_all(string = fasta$file_name, 
                                  pattern = "[\\(\\),\\[\\]]+") %>% 
  str_replace_all(pattern = " ", replacement = "_")

## FASTA clean sequences ----
## trim leading and trailing Ns
message("\t", "Cleaning FASTA sequences...")
message("\t\t", "Trimming leading/trailing Ns...")
fasta$Sequence <- str_remove_all(string = fasta$Sequence,
                                 pattern = "^N+|N+$")

## FASTA length < 400 ----
message("\t\t", "Removing FASTA sequences with LENGTH < 400 nt ...")
fasta$length <- str_length(string = fasta$Sequence)

idx <- which(fasta$length >= 400)
fasta <- fasta[idx,]

## FASTA count_N ----
message("\t\t", "Removing FASTA sequences with PERCENT_N > 50% ...")
fasta$count_N <- str_count(string = fasta$Sequence,
                           pattern = "N")

## FASTA percent_N ----
fasta$percent_N <- (fasta$count_N / fasta$length) * 100

## FASTA remove sequences with percent_N >= 50 ----
idx <- which(fasta$percent_N <= 50)
fasta <- fasta[idx,]

## ANALYSIS_INPUT ----
freq <- count(df = fasta, vars = "in_assembly")
freq <- rename_col(df = freq, 
                   old_name = "freq", 
                   new_name = "analysis_input")

## O_INPUT merge o_input and freq ----
col_names <- union(names(o_input), 
                   names(freq))

o_input <- merge(x = o_input, 
                 y = freq, 
                 by = "in_assembly", all = TRUE, sort = FALSE)[,col_names]

idx <- which(is.na(o_input$analysis_input))

o_input$analysis_input[idx] <- 0

## FASTA merge fasta and o_input ----
col_names <- union(names(fasta), 
                   names(o_input))

fasta <- merge(x = fasta, 
               y = o_input, 
               all.x = TRUE, sort = FALSE)[,col_names]

## FASTA analysis_input fas files ----
message("\n\t", "Generating ANALYSIS_INPUT files...")
Assembly_ID <- unique(fasta$Assembly_ID)

## OUT_FAS
out_fas <- vector()

## READ_1 & READ_2
Read_1 <- vector()
Read_2 <- vector()

## DIR_OUT
dir_out <- vector()

for (i in 1:length(Assembly_ID)) {
  idx <- which(fasta$Assembly_ID == Assembly_ID[i])
  sub_fasta <- fasta[idx,]
  
  dir_in <- paste0(analysis_in_dir, "/", i)
  dir_out[i] <- dir_in
  message("\t\t", "Assembly_ID: ", Assembly_ID[i])
  dir.create(dir_in, 
             showWarnings = FALSE)
  
  reads_dir <- paste0(dir_in, "/reads")
  dir.create(path = reads_dir, 
             showWarnings = FALSE)
  
  Header <- sub_fasta$seq.ID
  Sequence <- sub_fasta$Sequence
  
  fas <- build_fas(Header = Header, 
                   Sequence = Sequence)
  
  out_fas[i] <- paste0(dir_in, 
                       "/", Assembly_ID[i], ".fas")
  
  writeFasta(fdta = fas, out.file = out_fas[i])
  
  message("\t\t\t", out_fas[i])
  in_fq_1 <- unique(sub_fasta$Read_1)
  in_fq_2 <- unique(sub_fasta$Read_2)
  
  if (length(in_fq_1) > 0) {
    b_name <- basename(in_fq_1)
    b_name <- str_remove_all(string = b_name, pattern = "[\\(\\),\\[\\]]+") %>%
      str_replace_all(pattern = " ", replacement = "_")
    
    out_fq_1 <- paste0(reads_dir, "/", b_name)
    
    file.copy(from = in_fq_1, 
              to = out_fq_1)
    
    message(paste0("\t\t\t", out_fq_1, collapse = "\n"))
    
    Read_1[i] <- out_fq_1
    
  } else {
    Read_1[i] <- NA
  }
  
  if (length(in_fq_2) > 0) {
    b_name <- basename(in_fq_2)
    b_name <- str_remove_all(string = b_name, pattern = "[\\(\\),\\[\\]]+") %>%
      str_replace_all(pattern = " ", replacement = "_")
    
    out_fq_2 <- paste0(reads_dir, "/", b_name)
    
    file.copy(from = in_fq_2, 
              to = out_fq_2)
    
    message(paste0("\t\t\t", out_fq_2, collapse = "\n"))
    
    Read_2[i] <- out_fq_2
    
  } else {
    Read_2[i] <- NA
  }
}

## JOB_FILE ----
old_job_file <- grep(pattern = "job_file", 
                     x = file_list, value = TRUE)

if (length(old_job_file) > 0) {
  unlink(x = old_job_file, 
         recursive = TRUE, force = TRUE)
}

job_file <- data.frame(out_fas, Read_1, Read_2, dir_out)

colnames(job_file) <- c("IN_FASTA", "READ_1", "READ_2", "OUT_DIR")

date_time <- format(Sys.time(), "%Y-%m-%d_%I-%M_%p")

out_job_file <- paste0(analysis_in_dir, 
                       '/job_file_', date_time, ".txt")

message("\n\t", "Generating analysis input file list (job_file)...",
        "\n\t\t", out_job_file)

write.table(x = job_file, 
            file = out_job_file, 
            sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

## SAMPLE_INFO ----
sample_info <- o_input[col_nm_info]

sample_info <- count(df = sample_info)

sample_info$freq <- NULL

vars <- c("NGS_INSTRUMENT_ID", "NGS_RUN_ID", "NGS_METHOD_ID", "Assembly_ID", 
          "Assembly_N", "NGS_ISOLATE_SN", "NGS_ISOLATE",  "NGS_SN", "NGS_N", 
          "original_input", "analysis_input", "collection.date", "country", 
          "iso_alpha3_code", "host", "host.common", "type", "isolation.source", 
          "common_name.prefix","in_assembly", "Read_1", "Read_2")

# idx <- which(!colnames(o_input) %in% vars)
# colnames(o_input)[idx]

assembly_info <- o_input[vars]

## SAMPLE_INFO OUT save template file ----
message("\n\t","Saving SAMPLE_INFO_TABLE: ", in_info, "\n")

file <- in_info
sheet_name <- "Sample_Info"
x <- sample_info

wb <- loadWorkbook(file = file)

removeWorksheet(wb = wb, sheet = sheet_name)
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, 
               tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), 
             widths = "auto")

sheet_name <- "Assembly_Info"
x <- assembly_info

if (isTRUE(sheet_name %in% sheets(wb))) {
  removeWorksheet(wb = wb, sheet = sheet_name)
}

addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, 
               tableStyle = 'TableStyleMedium9')

setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), 
             widths = "auto")

sheet_name <- "Analysis_Input"
x <- job_file

if (isTRUE(sheet_name %in% sheets(wb))) {
  removeWorksheet(wb = wb, sheet = sheet_name)
}
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, 
               tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), 
             widths = "auto")

o_head <- fasta[c("seq.ID", "original_header")]
sheet_name <- "Original_Headers"
x <- o_head

if (isTRUE(sheet_name %in% sheets(wb))) {
  removeWorksheet(wb = wb, sheet = sheet_name)
}
addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, 
               tableStyle = 'TableStyleMedium9')
setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), 
             widths = "auto")

saveWorkbook(wb = wb, file = file, 
             overwrite = TRUE)

