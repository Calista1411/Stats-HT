# This function runs Bayesian estimation of Markov order on every sequence file 
# in a folder
# Arguments: 
  # folder: path to the folder containing sequence files
  # pattern: file format, defaults to text files
  # recursive: logical argument that controls whether to also search in 
             # sub-folders within the folder, defaults to F
  # alpha: Dirichlet prior hyperparameters, defaults to 0.5
  # k_max: highest Markov order to be evaluated, defaults to 5
  # prior_type: specifies type of prior, either "uniform" or "penalty"
  # save_csv: path and name of file to save the output into, defaults to NULL
# Output: 
  # A data frame where each row corresponds to one order for one input file, with
  # the following columns: model order, log-likelihood, log-evidence, log prior,
  # log posterior, posterior probability over orders, k*

run_bayes_folder <- function(
    folder,
    pattern   = "\\.txt$",
    recursive = FALSE,
    alpha     = 0.5,
    k_max     = 5,
    prior     = c("uniform", "penalty"),
    save_csv  = NULL
) {
  
  prior <- match.arg(prior)
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = recursive)

  results <- purrr::map_dfr(files, function(f) {
    out <- run_bayes_sequence(f, k_max = k_max, alpha = alpha,
                              prior_type = prior)
    tibble(
      file          = basename(f),
      k             = out$k,
      log_likelihood= out$log_likelihood,
      log_evidence  = out$log_evidence,
      log_prior     = out$log_prior,
      log_posterior = out$log_posterior,
      posterior     = out$posterior,
      k_star        = attr(out, "k_star")
    )
  })
  
  if (!is.null(save_csv)) readr::write_csv(results, save_csv)
  results
}