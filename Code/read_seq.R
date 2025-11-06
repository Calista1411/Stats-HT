# Function to read a DNA sequence from a file, remove invalid characters,
# and encodes each base as a numeric vector using the mapping: 
# A = 1, C = 2, G = 3, T = 4
# Arguments: 
  # file_path: path to the file containing the DNA sequence
# Output: 
  # A numeric vector representing the sequence 

read_seq <- function(file_path) {
  raw_seq <- readLines(file_path, warn = FALSE)
  raw_seq <- paste(raw_seq, collapse = "")
  sample <- unlist(strsplit(gsub("[^ACGT]", "", raw_seq), ""))
  as.numeric(chartr("ACGT", "1234", sample))
}