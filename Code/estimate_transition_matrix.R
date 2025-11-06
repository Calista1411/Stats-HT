# This function estimates the transition probability matrix of a discrete Markov 
# chain from an observed sequence
# Arguments:
  # seq: numeric vector representing the DNA sequence (A = 1, C = 2, G = 3, T = 4)
  # sample_space: number of possible states in the state space, defaults to 4 for DNA sequences
  # order: Markov order of the sequence to estimate
# Output: 
  # A matrix of dimension sample_space^order x sample_space, with elements 
  # representing the estimated transition probabilities

estimate_transition_matrix <- function(seq, sample_space = 4, order = 1) {
  n <- length(seq)
  rows <- sample_space^order
  counts <- matrix(0, nrow = rows, ncol = sample_space)
  
  for (i in (order+1):n) {
    prev_state <- 0
    for (j in 1:order) {
      prev_state <- prev_state * sample_space + (seq[i - order + j - 1] - 1)
    }
    prev_state <- prev_state + 1  
    counts[prev_state, seq[i]] <- counts[prev_state, seq[i]] + 1
  }
  
  # Normalise rows
  row_sums <- rowSums(counts)
  probs <- counts / ifelse(row_sums == 0, 1, row_sums)
  return(probs)
}