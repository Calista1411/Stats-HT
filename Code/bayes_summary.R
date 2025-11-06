# This function runs Bayesian estimation of Markov order across various 
# combinations of Dirichlet prior hyperparameters and prior types, and summarises
# the resulting k* for each sequence file 
# Arguments:
  # folder: path to the folder containing sequence files
  # alpha: Dirichlet prior hyperparameters, defaults to 0.5
  # priors: specifies types of prior, defaults to c("uniform" or "penalty")
  # k_max: highest Markov order to be evaluated, defaults to 5
  # save_csv: path and name of file to save the output into, defaults to NULL
# Output: 
  # A dataframe with columns representing the file name, alpha value used, 
  # prior type, k*, and the posterior probability corresponding to k*

bayes_summary <- function(
    folder,
    alphas   = c(0.1, 1, 2),
    priors   = c("uniform", "penalty"),
    k_max    = 5,
    save_csv = NULL
) {

  pars <- tidyr::crossing(alpha = alphas, prior = priors)
  
  # Run Bayesian estimation for each combination
  sens <- purrr::map_dfr(seq_len(nrow(pars)), function(i) {
    p <- pars[i, ]
    run_bayes_folder(folder,
                     alpha = p$alpha,
                     prior = p$prior,
                     k_max = k_max) %>%
      dplyr::mutate(alpha = p$alpha, prior = p$prior)
  })
  
  # Summarise k* for each sequence
  summary <- sens %>%
    dplyr::group_by(file, alpha, prior) %>%
    dplyr::summarise(
      k_MAP    = k[which.max(posterior)],
      post_MAP = max(posterior),
      .groups  = "drop"
    )
  
  if (!is.null(save_csv)) readr::write_csv(summary, save_csv)
  
  summary
}