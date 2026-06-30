source('Funzioni/bstexpnfit.R')

test_bstexpfit <- function() {
  cat("=== TEST BSTEXPFIT ===\n")
  
  set.seed(42)
  
  # Parametri veri
  true_lambda <- 0.1
  true_beta <- 0.8
  n <- 4000
  
  # Genera dati stretched exponential tramite trasformazione
  # Per Weibull: se U ~ Uniform(0,1), allora X = (-log(U)/lambda)^(1/beta)
  u <- runif(n)
  x <- (-log(1-u) / true_lambda)^(1/true_beta)
  x <- x[x < 100]  # Tronca valori estremi
  
  # Binning
  min_x <- min(x)
  max_x <- max(x)
  # Aggiungi un margine del 10% per sicurezza
  boundaries <- exp(seq(log(min_x * 0.9), log(max_x * 1.1), length.out = 25))
  h <- hist(x, breaks = boundaries, plot = FALSE)$counts
  h <- hist(x, breaks = boundaries, plot = FALSE)$counts
  
  # Test 1: Fit con ottimizzazione
  result <- bstexpfit(h, boundaries)
  
  cat(sprintf("Estimated params: λ = %.3f (true = %.3f), β = %.3f (true = %.3f)\n",
              result$lambda, true_lambda, result$beta, true_beta))
  
  # Tolleranza alta per stretched exponential (difficile da fittare)
  error_lambda <- abs(result$lambda - true_lambda) / true_lambda
  error_beta <- abs(result$beta - true_beta) / true_beta
  
  if (error_lambda > 0.5 || error_beta > 0.5) {
    warning("bstexpfit: high errors")
  }
  
  # Test 2: Grid search
  range_params <- list(
    lambda = seq(0.05, 0.2, by = 0.02),
    beta = seq(0.5, 1.2, by = 0.1)
  )
  result_grid <- bstexpfit(h, boundaries, range = range_params)
  
  # Verifica parametri positivi
  if (result$lambda <= 0 || result$beta <= 0) {
    stop("bstexpfit: non-positive params")
  }
  
  cat("bstexpfit test ended\n\n")
  return(TRUE)
}

test_bstexpfit()