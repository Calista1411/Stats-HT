# This function computes BICs and AICs and returns the Markov order with the 
# minimum AIC and BIC value
# Arguments: 
  # file_path: path to the text file containing the descendant sequence 
# Output: 
  # List containing the AIC and BIC values for each order, and the orders 
  # corresponding to minimum AIC and BIC

run_ic_tests <- function(file_path) {
  
  # Read and preprocess sequence (converts ACGT to 1234)
  raw_seq <- readLines(file_path, warn = FALSE)
  raw_seq <- paste(raw_seq, collapse = "")
  sample <- unlist(strsplit(gsub("[^ACGT]", "", raw_seq), ""))
  sample <- as.numeric(chartr("ACGT", "1234", sample))
  
  # Compute BIC values
  bic_scores <- c(
    BIC0 = BIC0(sample),
    BIC1 = BIC1(sample),
    BIC2 = BIC2(sample),
    BIC3 = BIC3(sample),
    BIC4 = BIC4(sample),
    BIC5 = BIC5(sample)
  )
  
  # Compute AIC values
  aic_scores <- c(
    AIC0 = AIC0(sample),
    AIC1 = AIC1(sample),
    AIC2 = AIC2(sample),
    AIC3 = AIC3(sample),
    AIC4 = AIC4(sample),
    AIC5 = AIC5(sample)
  )
  
  best_bic_order <- which.min(bic_scores) - 1
  best_aic_order <- which.min(aic_scores) - 1
  
  return(list(
    best_AIC_order = best_aic_order,
    best_BIC_order = best_bic_order,
    AICs = aic_scores,
    BICs = bic_scores
  ))
}

# Helper functions 

# BIC0: Zero-order (i.i.d.)
BIC0 <- function(sample) {
  counts <- table(factor(sample, levels = 1:4))
  proportions <- counts / length(sample)
  log_likelihood <- sum(counts * log(proportions))
  bic <- 3 * log(length(sample)) - 2 * log_likelihood
  return(bic)
}

# BIC1: First-order
BIC1 <- function(sample) {
  all_combinations <- generate_combinations(4, 2)
  combinations <- embed(sample, 2)[, 2:1]
  words <- apply(combinations, 1, paste, collapse = "")
  nab <- as.vector(table(factor(words, levels = all_combinations)))
  mat <- matrix(nab, nrow = 4, ncol = 4, byrow = TRUE)
  na_dot <- rep(rowSums(mat), rep(4, 4))
  valid <- nab > 0
  bic <- 12 * log(length(sample)) - 2 * sum(nab[valid] * log(nab[valid] / na_dot[valid]))
  return(bic)
}

# BIC2: Second-order
BIC2 <- function(sample) {
  all_combinations <- generate_combinations(4, 3)
  combinations <- embed(sample, 3)[, 3:1]
  words <- apply(combinations, 1, paste, collapse = "")
  nabc <- as.vector(table(factor(words, levels = all_combinations)))
  mat <- matrix(nabc, nrow = 16, ncol = 4, byrow = TRUE)
  nab_dot <- rep(rowSums(mat), rep(4, 16))
  valid <- nabc > 0
  bic <- 48 * log(length(sample)) - 2 * sum(nabc[valid] * log(nabc[valid] / nab_dot[valid]))
  return(bic)
}

# BIC3: Third-order
BIC3 <- function(sample) {
  all_combinations <- generate_combinations(4, 4)
  combinations <- embed(sample, 4)[, 4:1]
  words <- apply(combinations, 1, paste, collapse = "")
  nabcd <- as.vector(table(factor(words, levels = all_combinations)))
  if (length(nabcd) < 256) nabcd <- c(nabcd, rep(0, 256 - length(nabcd)))
  mat <- matrix(nabcd, nrow = 64, ncol = 4, byrow = TRUE)
  nabc_dot <- rep(rowSums(mat), rep(4, 64))
  valid <- nabcd > 0 & nabc_dot > 0
  bic <- 192 * log(length(sample)) - 2 * sum(nabcd[valid] * log(nabcd[valid] / nabc_dot[valid]))
  return(bic)
}

# BIC4: Fourth-order
BIC4 <- function(sample) {
  all_combinations <- generate_combinations(4, 5)
  combinations <- embed(sample, 5)[, 5:1]
  words <- apply(combinations, 1, paste, collapse = "")
  nabcde <- as.vector(table(factor(words, levels = all_combinations)))
  if (length(nabcde) < 1024) nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  mat <- matrix(nabcde, nrow = 256, ncol = 4, byrow = TRUE)
  nabcd_dot <- rep(rowSums(mat), rep(4, 256))
  valid <- nabcde > 0 & nabcd_dot > 0
  bic <- 768 * log(length(sample)) - 2 * sum(nabcde[valid] * log(nabcde[valid] / nabcd_dot[valid]))
  return(bic)
}

# BIC5: Fifth-order
BIC5 <- function(sample) {
  all_combinations <- generate_combinations(4, 6)
  combinations <- embed(sample, 6)[, 6:1]
  words <- apply(combinations, 1, paste, collapse = "")
  nabcdef <- as.vector(table(factor(words, levels = all_combinations)))
  if (length(nabcdef) < 4096) nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  mat <- matrix(nabcdef, nrow = 1024, ncol = 4, byrow = TRUE)
  nabcde_dot <- rep(rowSums(mat), rep(4, 1024))
  valid <- nabcdef > 0 & nabcde_dot > 0
  bic <- 3072 * log(length(sample)) - 2 * sum(nabcdef[valid] * log(nabcdef[valid] / nabcde_dot[valid]))
  return(bic)
}

# AIC0: Zero-order (i.i.d.)
AIC0<-function(sample){
  counts <- table(factor(sample, levels = 1:4))
  proportions <- counts / length(sample)
  log_likelihood <- sum(counts * log(proportions))
  aic <- 2*3 - 2 * log_likelihood
  return(aic)
}

# AIC1: 1st-order 
AIC1 <- function(sample) {
  all_combinations <- generate_combinations(4, 2)
  window_size <- 2
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  two_letter_count <- table(factor(words, levels = all_combinations))
  nab <- as.vector(two_letter_count)
  two_words_matrix <- matrix(nab, nrow = 4, ncol = 4,byrow=TRUE)
  na.=rep(rowSums(two_words_matrix),rep(4,4))
  valid_indices=nab>0
  nab_valid=nab[valid_indices]
  na._valid=na.[valid_indices]
  aic=2*12-2*sum(nab_valid*log(nab_valid/na._valid))
  return(aic)
}

# AIC2: 2nd-order
AIC2 <- function(sample) {
  all_combinations <- generate_combinations(4, 3)
  window_size <- 3
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  three_letter_count <- table(factor(words, levels = all_combinations))
  nabc <- as.vector(three_letter_count)
  three_words_matrix <- matrix(nabc, nrow = 4^2, ncol = 4,byrow=TRUE)
  nab.=rep(rowSums(three_words_matrix),rep(4,4^2))
  valid=nabc>0
  aic=2*48-2*sum(nabc[valid]*log(nabc[valid]/nab.[valid]))
  return(aic)
}

# AIC3: 3rd-order
AIC3 <- function(sample) {
  all_combinations <- generate_combinations(4, 4)
  window_size <- 4
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  four_letter_count <- table(factor(words, levels = all_combinations))
  nabcd <- as.vector(four_letter_count)
  if (length(nabcd) < 64) {
    nabcd <- c(nabcd, rep(0, 64 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  valid=nabcd>0 & nabc.>0
  aic=2*192-2*sum(nabcd[valid]*log(nabcd[valid]/nabc.[valid]))
  return(aic)
}

# AIC4: 4th-order 
AIC4 <-function(sample) {
  all_combinations <- generate_combinations(4, 5)
  window_size <- 5
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  five_letter_count <- table(factor(words, levels = all_combinations))
  nabcde <- as.vector(five_letter_count)
  if (length(nabcde) < 1024) { 
    nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  }
  five_words_matrix <- matrix(nabcde, nrow = 4^4, ncol = 4,byrow=TRUE)
  nabcd.=rep(rowSums(five_words_matrix),rep(4,4^4))
  valid=nabcde>0 & nabcd.>0
  aic=2*768-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(aic)
}

# AIC5: 5th-order
AIC5<- function(sample) {
  all_combinations <- generate_combinations(4, 6)
  window_size <- 6
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  six_letter_count <- table(factor(words, levels = all_combinations))
  nabcdef <- as.vector(six_letter_count)
  if (length(nabcdef) < 4096) { 
    nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  valid=nabcdef>0 & nabcde. >0
  aic=2*3072-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(aic)
}