source ("Funzioni/bexpnfit.R")

test_bexpnfit <- function() {
  cat("=== TEST BEXPNFIT ===\n")
  
  set.seed(111)
  
  # Genera dati da distribuzione esponenziale vera
  true_lambda <- 0.5
  true_bmin <- 1
  n <- 5000
  
  # Genera dati esponenziali troncati
  x <- rexp(n, rate = true_lambda) + true_bmin
  
  # Binning logaritmico
  boundaries <- exp(seq(log(0.5), log(max(x) + 1), length.out = 30))
  h <- hist(x, breaks = boundaries, plot = FALSE)$counts
  
  # Test 1: Fit con ottimizzazione
  result <- bexpnfit(h, boundaries, bmin = true_bmin)
  
  # Verifica accuratezza
  error <- abs(result$lambda - true_lambda) / true_lambda
  cat(sprintf("Lambda: estimate = %.3f, true = %.3f, error = %.1f%%\n",
              result$lambda, true_lambda, error * 100))
  
  if (error > 0.2) {
    warning("bexpnfit: error > 20%")
  }
  
  # Test 2: Fit con grid search
  lambda_range <- seq(0.1, 1.0, by = 0.05)
  result_grid <- bexpnfit(h, boundaries, lambda_range = lambda_range, bmin = true_bmin)
  
  # Verifica che grid search dia risultati simili
  if (abs(result$lambda - result_grid$lambda) / result$lambda > 0.1) {
    warning("bexpnfit: incoherence between optimization and grid search")
  }
  
  # Test 3: Log-likelihood deve essere finita e negativa
  if (!is.finite(result$loglik) || result$loglik >= 0) {
    stop("bexpnfit: log-likelihood not valid")
  }
  
  cat("bexpnfit test ended\n\n")
  return(TRUE)
}


test_bexpnfit()
