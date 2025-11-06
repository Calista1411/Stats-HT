# This function runs AIC and BIC on all files within a folder 
# Arguments:
  # folder_path: path to folder
  # pattern: file format
  # recursive: logical argument that controls whether to also search in 
             # sub-folders within the folder, defaults to F
  # save_csv_folder: path and name of file to save the output in, defaults to NULL
# Output: 
  # Dataframe containing the AIC and BIC values for each order and the order 
  # with minimum AIC and BIC for each sequence file in the folder

run_ic_folder <- function(
    folder_path,
    pattern   = "\\.txt$",
    recursive = FALSE,
    save_csv_folder  = NULL) {
  
  files <- list.files(folder_path, pattern = pattern, full.names = TRUE, recursive = recursive)
  
  rows <- list()
  
  for (f in files) {
    cat("\nProcessing:", basename(f), "\n")
    res <- tryCatch(
      run_ic_tests(f),
      error = function(e) { warning("Error in file ", basename(f), ": ", e$message); NULL }
    )
    if (is.null(res)) next
    
    row_list <- c(
      list(
        file = basename(f),
        best_AIC_order = res$best_AIC_order,
        best_BIC_order = res$best_BIC_order
      ),
      as.list(res$AICs),
      as.list(res$BICs)
    )
    rows[[length(rows) + 1]] <- as.data.frame(row_list, check.names = FALSE)
  }
  
  combined_df <- do.call(rbind, rows)
  rownames(combined_df) <- NULL
  
  if (!is.null(save_csv_folder)) {
    write.csv(combined_df, save_csv_folder, row.names = FALSE)
    message("Saved results to: ", normalizePath(save_csv_folder, winslash = "/"))
  }
  
  combined_df
}