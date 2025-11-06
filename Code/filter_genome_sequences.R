# This function filters and extracts coding DNA sequences with lengths that fall 
# in the specified ranges
# Arguments:
  # folder: path to folder containing CSV files with sequences to be filtered
  # windows: named list containing sequence length ranges of interest, corresponding 
           # to approximately 5,000 bases and 10,000 bases
  # pattern: file format, defaults to CSV
  # recursive: logical argument that controls whether to also search in 
             # sub-folders within the folder, defaults to F
  # output_dir: directory where filtered CSVs will be written, defaults to the input folder
# Output: None, CSV files containing the filtered sequences will be written and saved

filter_genome_sequences <- function(
    folder,
    windows   = list(`5000` = c(4980, 5020), `10000` = c(9700, 10300)),
    pattern   = "\\.csv$",
    recursive = FALSE,
    output_dir  = folder
) {
  stopifnot(dir.exists(folder))
  
  read_one <- function(f) {
    readr::read_csv(
      f,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ) %>%
      mutate(source_file = basename(f))
  }
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(files) == 0L) stop("No CSV files found in '", folder, "' (recursive = ", recursive, ").")
  
  all_data <- purrr::map_dfr(files, read_one) %>%
    dplyr::rename(
      Chr      = Chr,
      CDS_Name = `CDS Name`,
      Sequence = Sequence
    ) %>%
    dplyr::mutate(
      Sequence   = toupper(stringr::str_replace_all(Sequence, "[^ACGT]", "")),
      seq_length = nchar(Sequence)
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\nSummary of windows:\n")
  walk2(names(windows), windows, function(name, rng) {
    stopifnot(is.numeric(rng), length(rng) == 2, rng[1] <= rng[2])
    
    # Selects only the coding DNA sequences with lengths between the specified ranges
    df <- all_data %>%
      filter(!is.na(seq_length), seq_length >= rng[1], seq_length <= rng[2]) %>%
      select(Chr, CDS_Name, seq_length, source_file, Sequence)
    
    cat(sprintf("\nSequences around %s bp (%d–%d): %d sequences\n",
                name, rng[1], rng[2], nrow(df)))
    if (nrow(df) > 0) {
      print(df %>% select(Chr, CDS_Name, seq_length, source_file), n = min(10, nrow(df)))
      if (nrow(df) > 10) cat("(showing first 10)\n")
    }
    
    out_path <- file.path(output_dir, paste0("seqs_around_", name, ".csv"))
    write_csv(df, out_path)
    cat("→ Written:", out_path, "\n")
  })
  
  invisible(NULL)
}