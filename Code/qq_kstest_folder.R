# For each sequence file in a folder, this function:
# (i) estimates the null Markov model P with estimate_transition_matrix(), 
# (ii) simulates an empirical null distribution of LR statistics via 
     # generate_null_distribution(),
# (iii) compares the simulated statistics to the theoretical distribution using
      # a Q–Q plot and a K-S test,
# (iv) plots a histogram of the distribution of test statistics

# Arguments: 
  # folder: path to folder containing the sequence files
  # pattern: file format
  # recursive: logical argument that controls whether to also search in 
             # sub-folders within the folder, defaults to F
  # sample_space: number of possible states in the state space, defaults to 4 
                # for DNA sequences
  # null_order: Markov order under the null hypothesis of the LR test 
  # n: length of each simulated sequence, defaults to 10,000
  # nrep: number of simulated sequences to generate, defaults to 1,000
  # seed: random seed for reproducibility, defaults to 1 (used in this thesis)
  # save_csv: path and name of file to save the output into, defaults to NULL
# Output: 
  # A dataframe with one row per file, containing file name, number of df, mean
  # and variance of simulated statistics, theoretical moments and K-S test outputs 
  # Q-Q plots and histograms for each file

qq_kstest_folder <- function(
    folder,
    pattern      = "\\.txt$",
    recursive    = FALSE,
    sample_space = 4,
    null_order   = 1,         
    n            = 10000,      
    nrep         = 1000,      
    seed         = 1,
    save_csv      = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = recursive)
  
  k  <- null_order
  df <- 3 * (sample_space^(k + 1)) - 3 * (sample_space^k)
  message("Expected Chi-square(df) = ", df)
  
  rows <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    f <- files[i]
    base <- tools::file_path_sans_ext(basename(f))
    cat(sprintf("[%d/%d] %s\n", i, length(files), basename(f)))
    
    seq_int <- read_seq(f)
    if (!length(seq_int)) next
    
    # Estimate transition matrix under the null model
    P_null <- estimate_transition_matrix(seq_int, sample_space = sample_space, order = null_order)
    
    exp_rows <- sample_space^null_order
    if (nrow(P_null) != exp_rows) {
      warning("Transition matrix row mismatch for ", basename(f),
              " (got ", nrow(P_null), ", expected ", exp_rows, "). Skipping.")
      next
    }
    
    # Simulate empirical null distribution under the null model
    null_stats <- generate_null_distribution(
      P            = P_null,
      n            = n,
      nrep         = nrep,
      sample_space = sample_space,
      order        = null_order
    )
    
    # KS test: empirical vs theoretical chi-square(df)
    ks <- suppressWarnings(stats::ks.test(null_stats, "pchisq", df = df))
    
    op <- par(no.readonly = TRUE)
    
    # Q–Q plot
    par(family = "CMU Serif")
    stats::qqplot(
      stats::qchisq(ppoints(length(null_stats)), df = df),
      sort(null_stats),
      xlab = bquote("Theoretical Quantiles (" * chi^2 * "(" * .(df) * "))"),
      ylab = "Empirical Quantiles",
      main = "Q–Q Plot",
      font.main = 2,
      pch = 19, cex = 0.5, col = "red"
    )
    abline(0, 1, col = "purple", lwd = 2)
    
    # Histogram with theoretical density overlay
    par(family = "CMU Serif")
    hist(null_stats, breaks = 60, freq = FALSE, col = "#FFB6C1",
         border = "white", main = "Distribution of Test Statistics",
         font.main = 2, xlab = "Test Statistic")
    curve(dchisq(x, df = df), add = TRUE, col = "red", lwd = 2)
    legend("topright",
           legend = c("Empirical", bquote(chi^2 * "(" * .(df) * ")")),
           col = c("#FFB6C1", "red"), lwd = c(6, 2), bty = "n")
    
    par(op)
    
    rows[[i]] <- data.frame(
      file         = basename(f),
      df           = df,
      mean_emp     = mean(null_stats),
      var_emp      = var(null_stats),
      mean_chisq   = df,
      var_chisq    = 2 * df,
      ks_statistic = unname(ks$statistic),
      ks_p_value   = ks$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  res <- do.call(rbind, rows)
  if (!is.null(save_csv) && nrow(res)) readr::write_csv(res, save_csv)
  res
}