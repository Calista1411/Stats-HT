# This function calculates the number and percentage of base positions that 
# differ between the ancestral and descendant sequence
# Arguments: 
  # ancestral_seq: path to file containing the ancestral sequence
  # descendant_seq: path to file containing the descendant sequence
# Output:
  # List containing the number and percentage of base positions that differ 
  # between the ancestral and descendant sequence

compare_sequences <- function(ancestral_seq, descendant_seq) {
  
  ancestral_seq <- read_seq(ancestral_seq)
  descendant_seq <- read_seq(descendant_seq)
  
  # Ensure sequences are same length
  if (length(ancestral_seq) != length(descendant_seq)) {
    stop("Sequences have different lengths.")
  }
  
  # Count the number and percentage of changed positions
  num_changed <- sum(ancestral_seq != descendant_seq)
  pct_changed <- (num_changed / length(ancestral_seq)) * 100
  
  list(
    num_changed = num_changed,
    pct_changed = pct_changed
  )
}