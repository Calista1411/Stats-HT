# This function generates all possible word combinations of a given length from 
# a finite state space 
# Arguments: 
  # n: number of states in state space (4 for DNA sequences, corresponding to 
     # the four nucleotides)
  # length: number of letters in each word combination
# Output: 
  # A character vector containing all possible word combinations (n^length) of the 
  # specified length. Each element is a string representing one combination.
  # Example: generate all dinucleotide combinations for DNA with generate_combinations(4, 2)
              # [1] "11" "12" "13" "14" "21" "22" ... "44"

generate_combinations <- function(n, length) {
  seq <- 1:n
  result <- list()
  generate <- function(current, remaining) {
    if (remaining == 0) {
      result <<- c(result, list(paste(current, collapse = "")))
    } else {
      for (i in seq) {
        generate(c(current, i), remaining - 1)
      }
    }
  }
  generate(integer(0), length)
  unlist(result)
}