#!/usr/bin/Rscript

## PACKAGES  ----
required_pkgs <- c("plyr", "stringr", "microseq", "openxlsx", "tools", "gtools", "optparse")

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

message('RUNNING: 1_generate_sample_info_table.R',
        '\n\n\t', 'IN_DIR: ', in_dir)

## PROJECT_NAME ----
project_name <- basename(in_dir)

## FUNCTIONS ----
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
analysis_out_dir <- paste0(in_dir, "/3_analysis_output")
table2asn_in_dir <- paste0(in_dir, "/4_table2asn_input")
table2asn_out_dir <- paste0(in_dir, "/5_table2asn_output")

dir_list <- c(analysis_in_dir, 
              analysis_out_dir)

clear_and_create_dirs(dir_list = dir_list, 
                      create = TRUE)

dir_list <- c(table2asn_in_dir, 
              table2asn_out_dir)

clear_and_create_dirs(dir_list = dir_list, 
                      create = FALSE)

message("\n\t", "ANALYSIS_IN_DIR: ", analysis_in_dir)
message("\t", "ANALYSIS_OUT_DIR: ", analysis_out_dir)

## FILE_LIST ----
file_list <- list.files(path = in_dir, 
                        full.names = TRUE, recursive = TRUE)

## ORIGINAL_SEQUENCES o_file_list ----
o_file_list <- grep(pattern = "1_original_input", 
                    x = file_list, value = TRUE)

## ORIGINAL_SEQUENCES a_dir ----
a_dir <- dirname(o_file_list)
a_dir <- grep(pattern = "assembly$", 
                     x = a_dir, value = TRUE) %>% unique()

## O_INPUT ----
o_input <- data.frame()
for (i in 1:length(a_dir)) {
  m_dir <- dirname(path = a_dir[i])
  
  file_name <- list.files(path = m_dir, 
                          full.names = TRUE, recursive = TRUE)
  
  a_files <- data.frame(file_name)
  
  a_files$NGS_ISOLATE_SN <-  basename(file_name) %>% 
    str_extract(pattern = "^.+_S[0-9]+_") %>% 
    str_remove(pattern = '_$')
  
  NGS_ISOLATE_SN <- unique(a_files$NGS_ISOLATE_SN)
  
  df_out <- data.frame()
  
  for (j in 1:length(NGS_ISOLATE_SN)) {
    idx <- which(a_files$NGS_ISOLATE_SN == NGS_ISOLATE_SN[j])
    sub_a_files <- a_files[idx,]
    
    in_assembly <- grep_na(pattern = "assembly", 
                        x = sub_a_files$file_name, 
                        value = TRUE)
    
    Read_1 <- grep_na(pattern = ".+_R1_*.+$", 
                      x = sub_a_files$file_name, 
                      value = TRUE)
  
    Read_2 <- grep_na(pattern = ".+_R2_*.+$", 
                      x = sub_a_files$file_name, 
                      value = TRUE)
    
    df <- data.frame(in_assembly, Read_1, Read_2)
    
    df$m_dir <- m_dir
    
    df$a_dir <- a_dir[i]
    
    df$info <- str_remove(string = df$m_dir, 
                          pattern = "^.+/1_original_input/")
    
    split_parts <- str_split_fixed(string = df$info, 
                                   pattern = "/", n = 3) %>% data.frame()
    
    colnames(split_parts) <- c("NGS_INSTRUMENT_ID", 
                               "NGS_RUN_ID", 
                               "NGS_METHOD_ID")
    
    df <- cbind(df, split_parts)
    
    df$NGS_ISOLATE_SN <- NGS_ISOLATE_SN[j]
    
    df$NGS_ISOLATE <- str_remove(string = df$NGS_ISOLATE_SN, 
                                  pattern = "_S[0-9]+$")
    
    df$NGS_SN <- str_extract(string = df$NGS_ISOLATE_SN, 
                                        pattern = "_S[0-9]+$") %>% str_remove(pattern = "^_")
    
    df$NGS_N <- str_remove(string = df$NGS_SN, pattern = "^S")
    
    df_out <- rbind(df_out, df)
  }
  
  o_input <- rbind(o_input, df_out)
}

## O_INPUT sort ----
o_input <- o_input[mixedorder(x = o_input$NGS_METHOD_ID, 
                              decreasing = FALSE),]

o_input <- o_input[mixedorder(x = o_input$NGS_SN, 
                              decreasing = FALSE),]

o_input <- o_input[mixedorder(x = o_input$NGS_RUN_ID, 
                              decreasing = FALSE),]

rownames(o_input) <- NULL

## O_INPUT Assembly_N ----
message("\n\t", "Generating ASSEMBLY_IDs...")
max_len <- rownames(o_input) %>% str_length() %>% max()

o_input$Assembly_N <- str_pad(string = rownames(o_input), 
                              width = max_len, side = 'left', pad = '0')

o_input$Assembly_N <- paste0('Assembly_', o_input$Assembly_N)

## O_INPUT Assembly_ID ({Assembly_N}-{NGS_ISOLATE_SN}) ----
o_input$Assembly_ID <- paste0(o_input$Assembly_N, "-", 
                              o_input$NGS_ISOLATE_SN)

paste("\t\t", paste0(o_input$Assembly_N, ": "), 
      o_input$Assembly_ID, 
      o_input$NGS_RUN_ID, 
      o_input$NGS_METHOD_ID, "\n" , 
      sep = "\t") %>% message()

## O_INPUT out_o_input ----

o_input <- o_input %>% select(-c(m_dir, a_dir, info))


out_o_input <- paste0(analysis_in_dir, 
                      "/original_input-", project_name, ".txt" )

message("\t", "Writing OUT_O_INPUT table: ", out_o_input)

write.table(x = o_input, 
            file = out_o_input, 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

## SAMPLE_INFO ----
vars <- c("NGS_INSTRUMENT_ID", "NGS_RUN_ID", "NGS_ISOLATE_SN")

info <- plyr::count(df = o_input, vars = vars)
info$freq <- NULL

## SAMPLE_INFO collection.date ----
info$collection.date <- ""

## SAMPLE_INFO country ----
info$country <- ""

## SAMPLE_INFO Iso_alpha3_code ----
info$iso_alpha3_code <- ""

## SAMPLE_INFO Host ----
info$host <- ""

## SAMPLE_INFO Host.common ----
info$host.common <- ""

## SAMPLE_INFO Type ----
info$type <- ""

## SAMPLE_INFO isolation.source ----
info$isolation.source <- ""

## SAMPLE_INFO Common_name.prefix ----
info$common_name.prefix <- ""


## SAMPLE_INFO out_info ----
out_info <- paste0(analysis_in_dir, 
                   "/Sample_Info-", project_name, ".xlsx")

wb <- createWorkbook()
out_wb <- out_info
sheet_name <- 'Sample_Info'
x <- info

addWorksheet(wb = wb, sheetName = sheet_name)
writeDataTable(wb = wb, sheet = sheet_name, x = x, 
               colNames = TRUE, rowNames = FALSE, 
               tableStyle = 'TableStyleMedium9')

setColWidths(wb = wb, sheet = sheet_name, 
             cols = 1:ncol(x), widths = "auto")

saveWorkbook(wb = wb, file = out_wb, overwrite = TRUE)

message("\n\t","Writing SAMPLE_INFO_TABLE: ", out_wb, "\n")






