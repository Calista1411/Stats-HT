# This function runs sequential LR tests on all files in a folder 
# For each sequence file, this function 
# (i) reads the sequence with read_seq(),
# (ii) optionally builds empirical null distributions for the Order 1 vs Order 2 
     # and Order 2 vs Order 3 tests using generate_null_distribution()
# (iii) and runs LRTs via run_lrt_by_sections() on the full sequence or 
      # user-specified sections
# Arguments: 
  # folder: path to folder containing input sequence files
  # pattern: file format, defaults to text files 
  # recursive: logical argument that controls whether to also search in 
             # sub-folders within the folder, defaults to F
  # alpha: significance level, defaults to 0.05 
  # n_sections: optional input: if provided, splits each sequence into (equal) 
              # sections for separate LRTs
  # section_lengths: optional input: integer vector of section lengths, 
                   # overrides n_sections if provided
  # use_empirical: logical argument to construct empirical null distributions, 
                 # defaults to F
  # null_n: length of each simulated chain when building empirical nulls, 
          # defaults to 10,000
  # null_nrep_p12: number of replicates for the Order 1 null, defaults to 1,000
  # null_nrep_p23: number of replicates for the Order 2 null, defaults to 1,000
  # save_csv: path and name of file to save the output into, defaults to NULL
# Output
  # A data frame with one row per file × section including estimated order and 
  # LR test p-values 

run_lrt_folder <- function(
    folder,
    pattern = "\\.txt$",
    recursive = FALSE,
    alpha = 0.05,
    n_sections = NULL,
    section_lengths = NULL,
    use_empirical = FALSE,    
    null_n = 10000,            
    null_nrep_p12 = 1000,      
    null_nrep_p23 = 1000,       
    save_csv = NULL
) {
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (!length(files)) stop("No sequence files found in folder.")
  
  df_list <- list()
  
  for (f in files) {
    cat("\nProcessing:", basename(f), "\n")
    
    seq_num <- read_seq(f)
    if (use_empirical) {
      P1 <- estimate_transition_matrix(seq_num, sample_space = 4, order = 1)
      P2 <- estimate_transition_matrix(seq_num, sample_space = 4, order = 2)
    }
    
    if (use_empirical) {
      null12 <- generate_null_distribution(P1, n = null_n, nrep = null_nrep_p12,
                                           sample_space = 4, order = 1)
      null23 <- generate_null_distribution(P2, n = null_n, nrep = null_nrep_p23,
                                           sample_space = 4, order = 2)
      emp_null <- list(p12 = null12, p23 = null23)
    } else {
      emp_null <- list(p12 = NULL, p23 = NULL)
    }
    
    res <- tryCatch(
      run_lrt_by_sections(
        file_path = f,
        alpha = alpha,
        n_sections = n_sections,
        section_lengths = section_lengths,
        emp_null = emp_null
      ),
      error = function(e) { warning("Error in file ", basename(f), ": ", e$message); NULL }
    )
    if (is.null(res)) next
    
    get_safe <- function(x, name) if (!is.null(x) && name %in% names(x)) x[[name]] else NA
    add_row <- function(sec_name, r) {
      pv <- r$pvalues
      data.frame(
        file = basename(f),
        section = sec_name,
        order_estimate = r$order,
        p01        = get_safe(pv, "p01"),
        p12_chisq  = get_safe(pv, "p12.p_chisq"),
        p12_emp    = get_safe(pv, "p12.p_emp"),
        p23_chisq  = get_safe(pv, "p23.p_chisq"),
        p23_emp    = get_safe(pv, "p23.p_emp"),
        p34        = get_safe(pv, "p34"),
        p45        = get_safe(pv, "p45"),
        stringsAsFactors = FALSE
      )
    }
    
    if (is.list(res) && all(grepl("^Section_", names(res)))) {
      for (sn in names(res)) df_list[[length(df_list)+1]] <- add_row(sn, res[[sn]])
    } else {
      df_list[[length(df_list)+1]] <- add_row("Full sequence", res)
    }
  }
  
  combined_df <- do.call(rbind, df_list)
  rownames(combined_df) <- NULL
  
  if (!is.null(save_csv)) {
    readr::write_csv(combined_df, save_csv)
    message("Saved results to: ", normalizePath(save_csv, winslash = "/"))
  }
  
  combined_df
}