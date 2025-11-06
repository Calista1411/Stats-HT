# This function runs sequential LRT tests to determine Markov order
# Arguments: 
  # sample: numeric vector representing the sequence (A = 1, C = 2, G = 3, T = 4) 
  # alpha: significance level, defaults to 0.05
  # emp_null: list of vectors of empirical null test statistics for LR tests of 
            # Order 1 vs Order 2 and Order 2 vs Order 3
# Output
  # List containing the p-values for each LR test and the estimated order

run_lr_tests <- function(sample, alpha = 0.05, emp_null = list(p12 = NULL, p23 = NULL)) {
  
  # Order 0 vs Order 1
  p01 <- tryCatch(lrt_0_vs_1(sample), error = function(e) NA_real_)
  if (is.na(p01) || p01 >= alpha) {
    return(list(order = 0, pvalues = c(p01 = p01)))
  }
  
  # Order 1 vs Order 2
  p12_res <- tryCatch(lrt_1_vs_2(sample, empirical_null_dist = emp_null$p12),
                      error = function(e) list(p_chisq = NA_real_, p_emp = NA_real_))
  if (is.na(p12_res$p_chisq) || p12_res$p_chisq >= alpha) {
    return(list(order = 1,
                pvalues = c(p01 = p01,
                            "p12.p_chisq" = p12_res$p_chisq,
                            "p12.p_emp"   = p12_res$p_emp)))
  }
  
  # Order 2 vs Order 3
  p23_res <- tryCatch(lrt_2_vs_3(sample, empirical_null_dist = emp_null$p23),
                      error = function(e) list(p_chisq = NA_real_, p_emp = NA_real_))
  if (is.na(p23_res$p_chisq) || p23_res$p_chisq >= alpha) {
    return(list(order = 2,
                pvalues = c(p01 = p01,
                            "p12.p_chisq" = p12_res$p_chisq,
                            "p12.p_emp"   = p12_res$p_emp,
                            "p23.p_chisq" = p23_res$p_chisq,
                            "p23.p_emp"   = p23_res$p_emp)))
  }
  
  # Order 3 vs Order 4
  p34 <- tryCatch(lrt_3_vs_4(sample), error = function(e) NA_real_)
  if (is.na(p34) || p34 >= alpha) {
    return(list(order = 3,
                pvalues = c(p01 = p01,
                            "p12.p_chisq" = p12_res$p_chisq,
                            "p12.p_emp"   = p12_res$p_emp,
                            "p23.p_chisq" = p23_res$p_chisq,
                            "p23.p_emp"   = p23_res$p_emp,
                            p34 = p34)))
  }
  
  # Order 4 vs Order 5
  p45 <- tryCatch(lrt_4_vs_5(sample), error = function(e) NA_real_)
  if (is.na(p45) || p45 >= alpha) {
    return(list(order = 4,
                pvalues = c(p01 = p01,
                            "p12.p_chisq" = p12_res$p_chisq,
                            "p12.p_emp"   = p12_res$p_emp,
                            "p23.p_chisq" = p23_res$p_chisq,
                            "p23.p_emp"   = p23_res$p_emp,
                            p34 = p34, p45 = p45)))
  }
  
  return(list(order = 5,
              pvalues = c(p01 = p01,
                          "p12.p_chisq" = p12_res$p_chisq,
                          "p12.p_emp"   = p12_res$p_emp,
                          "p23.p_chisq" = p23_res$p_chisq,
                          "p23.p_emp"   = p23_res$p_emp,
                          p34 = p34, p45 = p45)))
}

# Helper functions

# LR test (Order 0 vs Order 1)
lrt_0_vs_1 <- function(sample) {
  words <- apply(embed(sample, 2)[, 2:1], 1, paste, collapse = "")
  freq2 <- table(factor(words, levels = generate_combinations(4, 2)))
  nab <- as.vector(freq2)
  matrix2 <- matrix(nab, nrow = 4, byrow = TRUE)
  na. <- rep(rowSums(matrix2), each = 4)
  counts <- tabulate(sample, nbins = 4)
  nb <- rep(counts, 4)
  n. <- rep(sum(counts), 16)
  expected <- (na. * nb) / n.
  valid <- !is.na(expected) & expected > 0
  test <- sum((nab[valid] - expected[valid])^2 / expected[valid])
  pchisq(test, df = 9, lower.tail = FALSE)
}

# LR test (Order 1 vs Order 2)
lrt_1_vs_2 <- function(sample, empirical_null_dist = NULL) {
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
  expected_check <- expected >= 5 
  print(expected_check)
  valid <- !is.na(expected) & expected > 0
  test <- sum((nabc[valid] - expected[valid])^2 / expected[valid])
  # Empirical p-value
  p_emp <- if (!is.null(empirical_null_dist)) {
    mean(empirical_null_dist >= test)
  } else {
    NA
  }
  p_chisq <- pchisq(test, df = 36, lower.tail = FALSE)
  return(list(p_chisq = p_chisq, p_emp = p_emp))
}

# LR test (Order 2 vs Order 3)
lrt_2_vs_3 <- function(sample, empirical_null_dist = NULL) {
  words4 <- apply(embed(sample, 4)[, 4:1], 1, paste, collapse = "")
  freq4 <- table(factor(words4, levels = generate_combinations(4, 4)))
  nabcd <- as.vector(freq4)
  matrix4 <- matrix(nabcd, nrow = 64, byrow = TRUE)
  nabc. <- rep(rowSums(matrix4), each = 4)
  words3 <- apply(embed(sample[-1], 3)[, 3:1], 1, paste, collapse = "")
  freq3 <- table(factor(words3, levels = generate_combinations(4, 3)))
  nbc. <- rep(rep(rowSums(matrix(freq3, nrow = 16, byrow = TRUE)), each = 4), 1)
  nbcd <- rep(as.vector(freq3), 4)
  expected <- (nabc. * nbcd) / nbc.
  expected_check <- expected >= 5
  print(expected_check)
  valid <- !is.na(expected) & expected > 0
  test <- sum((nabcd[valid] - expected[valid])^2 / expected[valid])
  # Empirical p-value
  p_emp <- if (!is.null(empirical_null_dist)) {
    mean(empirical_null_dist >= test)
  } else {
    NA
  }
  p_chisq <- pchisq(test, df = 144, lower.tail = FALSE)
  return(list(p_chisq = p_chisq, p_emp = p_emp))
}

# LR test (Order 3 vs order 4)
lrt_3_vs_4 <- function(sample) {
  words5 <- apply(embed(sample, 5)[, 5:1], 1, paste, collapse = "")
  freq5 <- table(factor(words5, levels = generate_combinations(4, 5)))
  nabcde <- as.vector(freq5)
  matrix5 <- matrix(nabcde, nrow = 256, byrow = TRUE)
  nabcd. <- rep(rowSums(matrix5), each = 4)
  words4 <- apply(embed(sample[-1], 4)[, 4:1], 1, paste, collapse = "")
  freq4 <- table(factor(words4, levels = generate_combinations(4, 4)))
  nbcd. <- rep(rep(rowSums(matrix(freq4, nrow = 64, byrow = TRUE)), each = 4), 1)
  nbcde <- rep(as.vector(freq4), 4)
  expected <- (nabcd. * nbcde) / nbcd.
  valid <- !is.na(expected) & expected > 0
  test <- sum((nabcde[valid] - expected[valid])^2 / expected[valid])
  pchisq(test, df = 576, lower.tail = FALSE)
}

# LR test (Order 4 vs Order 5)
lrt_4_vs_5 <- function(sample) {
  words6 <- apply(embed(sample, 6)[, 6:1], 1, paste, collapse = "")
  freq6 <- table(factor(words6, levels = generate_combinations(4, 6)))
  nabcdef <- as.vector(freq6)
  matrix6 <- matrix(nabcdef, nrow = 1024, byrow = TRUE)
  nabcde. <- rep(rowSums(matrix6), each = 4)
  words5 <- apply(embed(sample[-1], 5)[, 5:1], 1, paste, collapse = "")
  freq5 <- table(factor(words5, levels = generate_combinations(4, 5)))
  nbcde. <- rep(rep(rowSums(matrix(freq5, nrow = 256, byrow = TRUE)), each = 4), 1)
  nbcdef <- rep(as.vector(freq5), 4)
  expected <- (nabcde. * nbcdef) / nbcde.
  valid <- !is.na(expected) & expected > 0
  test <- sum((nabcdef[valid] - expected[valid])^2 / expected[valid])
  pchisq(test, df = 2304, lower.tail = FALSE)
}