# This function runs Bayesian estimation of Markov order on a sequence file, by 
# fitting Markov chain models of various orders and computing, for each order, 
# the Dirichlet-multinomial model evidence and the resulting posterior over orders
# Arguments: 
  # file_path: path to text file containing the DNA sequence
  # k_max: highest Markov order to evaluate, defaults to 5
  # alpha: Dirichlet prior hyperparameters, defaults to 0.5
  # prior_type: takes on either "uniform" or "penalty"
  # return_counts: logical vector to return transition count matrices, defaults to F
# Output: 
  # A dataframe with one row per order and the following columns: model order, 
  # log-likelihood at MLE, log-evidence, log prior, log posterior, posterior 
  # probability over orders, and the chosen order (k*)

run_bayes_sequence <- function(
    file_path,
    k_max         = 5,
    alpha         = 0.5,                
    prior_type    = c("uniform", "penalty"),
    return_counts = FALSE) {
  
  prior_type <- match.arg(prior_type)
  x <- read_seq(file_path)
  
  res <- lapply(0:k_max, function(k) {
    counts    <- count_ngram_transitions_int(x, k)
    log_ev    <- log_evidence_dirichlet_multinomial(counts, alpha)
    log_like  <- log_likelihood_mle(counts)
    log_prior <- log_model_prior(k, k_max, prior_type)
    list(k=k, counts=if (return_counts) counts else NULL,
         log_evidence=log_ev, log_likelihood=log_like, log_prior=log_prior)
  })
  
  # Posterior over orders (log(Pr(Mk|D)))
  log_post_unnorm <- vapply(res, function(z) z$log_evidence + z$log_prior, numeric(1))
  m <- max(log_post_unnorm)
  logZ <- m + log(sum(exp(log_post_unnorm - m)))
  post <- exp(log_post_unnorm - logZ)
  
  out <- data.frame(
    k             = 0:k_max,
    log_likelihood= vapply(res, `[[`, numeric(1), "log_likelihood"),
    log_evidence  = vapply(res, `[[`, numeric(1), "log_evidence"),
    log_prior     = vapply(res, `[[`, numeric(1), "log_prior"),
    log_posterior = log_post_unnorm - logZ,
    posterior     = post
  )
  attr(out, "k_star") <- out$k[which.max(out$posterior)]
  if (return_counts) attr(out, "counts") <- lapply(res, `[[`, "counts")
  out
}

# Helper functions

# Enumerate all contexts of length k (to label rows of transition matrix)
# For order-k Markov chain, there are 4^k contexts 
contexts_k <- function(k) {
  if (k == 0) return("")
  ctx_num <- generate_combinations(4, length = k)
  # map 1,2,3,4 -> A,C,G,T 
  chartr("1234", "ACGT", ctx_num)
}

# Count transitions n_{context -> next state}; returns matrix with 
# rows = contexts, cols (next state) = A/C/G/T
# x: input sequence as an integer vector in {1,2,3,4}
count_ngram_transitions_int <- function(x, k) {
  x <- as.integer(x)
  n <- length(x)
  contexts <- contexts_k(k) # row labels 
  n_ctx <- max(1, 4^k) 
  counts <- matrix(0, nrow = n_ctx, ncol = 4,
                   dimnames = list(contexts, c("A","C","G","T"))) 
  if (n == 0 || n <= k) return(counts)
  
  if (k == 0) {
    # For order 0, count no. of A/C/G/T directly
    tab <- tabulate(x[!is.na(x)], nbins = 4)
    counts[1, ] <- tab
    return(counts)
  }
  
  base4pow <- 4^(k - 1)
  code <- NA_integer_
  
  if (all(!is.na(x[1:k]))) {
    code <- 0
    for (r in 1:k) code <- code * 4 + (x[r] - 1)
  }
  
  for (t in (k + 1):n) {
    if (is.na(x[t - 1])) {
      code <- NA_integer_
    } else if (is.na(code)) {
      window <- x[(t - k):(t - 1)]
      if (all(!is.na(window))) {
        code <- 0
        for (r in seq_len(k)) code <- code * 4 + (window[r] - 1)
      }
    } else {
      leftmost <- x[t - k]
      if (is.na(leftmost)) {
        code <- NA_integer_
      } else {
        code <- (code %% base4pow) * 4 + (x[t - 1] - 1)
      }
    }

    if (!is.na(code) && !is.na(x[t])) {
      counts[code + 1, x[t]] <- counts[code + 1, x[t]] + 1
    }
  }
  counts
}

# Log model evidence for Dirichlet-multinomial per context row, 
# summed across rows (log(Pr(D|Mk)))
# Counts: matrix (4^k x 4), alpha: prior hyperparameters, either scalar, 
# length-4 vector, or matrix-like counts
log_evidence_dirichlet_multinomial <- function(counts, alpha = 0.5) {
  if (length(alpha) == 1) {
    alpha <- rep(alpha, 4)
  }
  if (is.vector(alpha) && length(alpha) == 4) {
    alpha <- matrix(rep(alpha, each = nrow(counts)), nrow = nrow(counts), byrow = FALSE)
  }
  stopifnot(all(dim(alpha) == dim(counts)))
  
  row_sum <- function(a_row, n_row) {
    lgamma(sum(a_row)) - sum(lgamma(a_row)) + sum(lgamma(a_row + n_row)) - lgamma(sum(a_row + n_row))
  }
  total <- 0
  for (r in seq_len(nrow(counts))) total <- total + row_sum(alpha[r, ], counts[r, ]) # add up for all rows
  total
}

# Number of free parameters for order k: 3 * 4^k 
num_params <- function(k) 3 * (4^k)

# Log prior over model orders: type = "uniform" or "penalty" (log(Pr(Mk)))
log_model_prior <- function(k, k_max, type = c("uniform", "penalty")) {
  type <- match.arg(type)
  if (type == "uniform") return(-log(k_max + 1))
  logw <- vapply(0:k_max, function(j) -1 * num_params(j), numeric(1))
  logw_k <- logw[k + 1]
  # Normalise with log-sum-exp
  m <- max(logw)
  lse <- m + log(sum(exp(logw - m)))
  logw_k - lse
}

# Log-likelihood at MLE
log_likelihood_mle <- function(counts) {
  total <- 0
  for (r in seq_len(nrow(counts))) {
    n_row <- counts[r, ]
    n_tot <- sum(n_row)
    if (n_tot == 0) next
    p_hat <- n_row / n_tot
    total <- total + sum(ifelse(n_row > 0, n_row * log(p_hat), 0))
  }
  total
}