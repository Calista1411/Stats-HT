# This function generates a random ancestral DNA sequence of length n
# Each nucleotide (A, C, G, T) is drawn independently with the specified probabilities
# Arguments:
  # n: length of sequence to generate
  # probs: Numeric vector of nucleotide probabilities. Defaults to equal probabilities
  # seed: Random seed for reproducibility, defaults to 1
# Output:
  # A single character string representing the generated DNA sequence

generate_ancestral_seq <- function(
    n,
    probs = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    seed  = 1
) {
  set.seed(seed)
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size = n, replace = TRUE, prob = probs[bases])
  paste(sequence, collapse = "")
}