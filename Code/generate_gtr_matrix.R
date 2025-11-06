# This function generates a normalised GTR rate matrix 
# Arguments: 
  # a,b,c,d,e,f: exchangeability parameters corresponding to the six pairs of nucleotides;  
                # defaults are the values used in this thesis: a = 1.0, b = 2.0, 
                # c = 1.5, d = 1.2, e = 0.9, f = 0.8
  # pi_A,pi_C,pi_G,pi_T: stationary base frequencies for A, C, G and T; 
                       # defaults are the values used in this thesis: all 0.25
# Output: A 4 by 4 matrix representing the normalised GTR rate matrix with rows 
        # and columns labelled by nucleotides ("A", "C", "G", "T")

generate_gtr_matrix <- function(a = 1.0, b = 2.0, c = 1.5,
                                d = 1.2, e = 0.9, f = 0.8,
                                pi_A = 0.25, pi_C = 0.25, 
                                pi_G = 0.25, pi_T = 0.25) {
  
  # Stationary distribution vector
  pi <- c(pi_A, pi_C, pi_G, pi_T)
  names(pi) <- c("A", "C", "G", "T")
  
  # Create unnormalised Q matrix
  Q <- matrix(0, nrow = 4, ncol = 4)
  rownames(Q) <- colnames(Q) <- c("A", "C", "G", "T")
  
  # Off-diagonal entries
  Q["A", "C"] <- a * pi["C"]
  Q["A", "G"] <- b * pi["G"]
  Q["A", "T"] <- c * pi["T"]
  
  Q["C", "A"] <- a * pi["A"]
  Q["C", "G"] <- d * pi["G"]
  Q["C", "T"] <- e * pi["T"]
  
  Q["G", "A"] <- b * pi["A"]
  Q["G", "C"] <- d * pi["C"]
  Q["G", "T"] <- f * pi["T"]
  
  Q["T", "A"] <- c * pi["A"]
  Q["T", "C"] <- e * pi["C"]
  Q["T", "G"] <- f * pi["G"]
  
  # Diagonal entries are set such that each row sums to zero
  for (i in 1:4) {
    Q[i, i] <- -sum(Q[i, -i])
  }
  
  # Normalise to ensure expected rate of substitution per unit time is 1
  expected_rate <- -sum(pi * diag(Q))
  Q_normalised <- Q / expected_rate
  
  return(Q_normalised)
}