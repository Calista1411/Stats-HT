# This function simulates replicate Markov chains from a specified transition
# matrix to approximate the empirical null distribution of the LR test statistic
# Arguments: 
  # P: transition probability matrix for the Markov chain, from estimate_transition_matrix()
  # n: length of each simulated sequence, defaults to 10,000
  # nrep: number of simulated sequences to generate, defaults to 1,000
  # sample_space: number of possible states in the state space, defaults to 4 for 
                # DNA sequences
  # order: Markov order of the sequence to estimate
# Output: 
  # A numeric vector of length n, containing the simulated LR test statistics

generate_null_distribution <- function(P, n = 10000, nrep = 1000, sample_space = 4, order) {
  results <- numeric(nrep)
  for (i in seq_len(nrep)) {
    results[i] <- simulation(P, n, sample_space, order)
  }
  return(results)
}

# Helper functions 

# Generates a random Markov chain from the specified transition matrix P from 
# estimate_transition_matrix() 
simulate_chain <- function(P, n, sample_space = 4, order) {
  chain <- numeric(n)
  chain[1:order] <- sample(1:sample_space, order, replace = TRUE)
  
  for (i in (order+1):n) {
    prev_state <- 0
    for (j in 1:order) {
      prev_state <- prev_state * sample_space + (chain[i - order + j - 1] - 1)
    }
    prev_state <- prev_state + 1  
    
    chain[i] <- sample(1:sample_space, 1, prob = P[prev_state, ])
  }
  return(chain)
}

# Computes the LR test statistic for a given Markov chain
simulation <- function(P, n = 10000, sample_space = 4, order) {
  sample <- simulate_chain(P, n, sample_space, order)
  
  if (order == 2) {
    words4 <- apply(embed(sample, 4)[, 4:1], 1, paste, collapse = "")
    freq4 <- table(factor(words4, levels = generate_combinations(sample_space, 4)))
    nabcd <- as.vector(freq4)
    matrix4 <- matrix(nabcd, nrow = sample_space^3, byrow = TRUE)
    nabc. <- rep(rowSums(matrix4), each = sample_space)
    
    words3 <- apply(embed(sample[-1], 3)[, 3:1], 1, paste, collapse = "")
    freq3 <- table(factor(words3, levels = generate_combinations(sample_space, 3)))
    nbc. <- rep(rep(rowSums(matrix(freq3, nrow = sample_space^2, byrow = TRUE)), each = sample_space), 1)
    nbcd <- rep(as.vector(freq3), sample_space)
    
    expected <- (nabc. * nbcd) / nbc.
    valid <- !is.na(expected) & expected > 0
    test <- sum((nabcd[valid] - expected[valid])^2 / expected[valid])
    
  } else if (order == 1) {
    
    words3 <- apply(embed(sample, 3)[, 3:1], 1, paste, collapse = "")
    freq3 <- table(factor(words3, levels = generate_combinations(4, 3)))
    nabc <- as.vector(freq3)
    matrix3 <- matrix(nabc, nrow = 16, byrow = TRUE)
    nab. <- rep(rowSums(matrix3), each = 4)
    words2 <- apply(embed(sample[-1], 2)[, 2:1], 1, paste, collapse = "")
    freq2 <- table(factor(words2, levels = generate_combinations(4, 2)))
    nb. <- rep(rep(rowSums(matrix(freq2, nrow = 4, byrow = TRUE)), each = 4), 1)
    nbc <- rep(as.vector(freq2), 4)
    expected <- (nab. * nbc) / nb.
    valid <- !is.na(expected) & expected > 0
    test <- sum((nabc[valid] - expected[valid])^2 / expected[valid])
    
  } else {
    stop("Only order = 1 or 2 currently supported.")
  }
  
  return(test)
}