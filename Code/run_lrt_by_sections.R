# This function run sequential LR tests to estimate Markov order 
# either for the entire sequence or for predefined or equal-length sections
# Arguments: 
  # file_path: path to text file containing descendant sequence
  # alpha: significance level, defaults to 0.05
  # n_sections: optional input: if provided, splits each sequence into (equal) 
              # sections for separate LRTs
  # section_lengths: optional input: integer vector of section lengths, 
                   # overrides n_sections if provided
  # emp_null: list of vectors of empirical null test statistics for LR tests of 
            # Order 1 vs Order 2 and Order 2 vs Order 
# Output: 
  # A list of results, where each element corresponds to a section (or the entire
  # sequence if not split), containing the estimated Markov order and the p-values

run_lrt_by_sections <- function(file_path, alpha = 0.05, 
                                n_sections = NULL, 
                                section_lengths = NULL, 
                                emp_null = list(p12 = NULL, p23 = NULL)) {
  
  sample <- read_seq(file_path)
  seq_len <- length(sample)
  results <- list()
  
  if (!is.null(section_lengths)) {
    if (sum(section_lengths) > seq_len) {
      stop("Sum of section_lengths is longer than the sequence length.")
    }
    start <- 1L
    for (i in seq_along(section_lengths)) {
      end   <- min(start + section_lengths[i] - 1L, seq_len)
      chunk <- sample[start:end]
      results[[i]] <- run_lr_tests(chunk, alpha = alpha, emp_null = emp_null)
      start <- end + 1L
    }
    names(results) <- paste0("Section_", seq_along(section_lengths))
    
  } else {
    if (is.null(n_sections)) n_sections <- 1L
    chunk_size <- ceiling(seq_len / n_sections)
    for (i in seq_len(n_sections)) {
      start <- (i - 1L) * chunk_size + 1L
      end   <- min(i * chunk_size, seq_len)
      chunk <- sample[start:end]
      results[[i]] <- run_lr_tests(chunk, alpha = alpha, emp_null = emp_null)
    }
    names(results) <- paste0("Section_", seq_along(results))
  }
  
  if (length(results) == 1) return(results[[1]])
  results
}