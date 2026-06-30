# TEST PER BPLCUTFIT
# Power law con cutoff esponenziale

source ("Funzioni/bplcutfit.R")

test_bplcutfit <- function() {
  cat("=== TEST BPLCUTFIT ===\n")
  
  set.seed(42)
  
  # Parametri veri
  true_alpha <- 2.5
  true_lambda <- 0.01
  true_bmin <- 1
  n <- 5000
  
  # Genera dati power law con cutoff 
  x <- numeric(0)
  while(length(x) < n) {
    # Genera da power law
    u <- runif(1)
    candidate <- true_bmin * (1 - u)^(-1/(true_alpha - 1))
    # Accetta con probabilità exp(-lambda * candidate)
    if(runif(1) < exp(-true_lambda * candidate)) {
      x <- c(x, candidate)
    }
  }
  x <- x[1:n]
  
  # Binning
  boundaries <- exp(seq(log(0.5), log(max(x) + 1), length.out = 25))
  h <- hist(x, breaks = boundaries, plot = FALSE)$counts
  
  # Test 1: Fit con ottimizzazione
  result <- bplcutfit(h, boundaries, bmin = true_bmin)
  
  cat(sprintf("Alpha: estimate = %.3f, true = %.3f\n", 
              result$alpha, true_alpha))
  cat(sprintf("Lambda: estimate = %.4f, true = %.4f\n",
              result$lambda, true_lambda))
  
  # Test 2: Grid search
  range_alpha <- seq(2.0, 3.0, by = 0.2)
  range_lambda <- seq(0.005, 0.02, by = 0.003)
  result_grid <- bplcutfit(h, boundaries, 
                           range_alpha = range_alpha,
                           range_lambda = range_lambda,
                           bmin = true_bmin)
  
  # Verifica constraints
  if (result$alpha <= 1) {
    stop("bplcutfit: alpha <= 1")
  }
  if (result$lambda <= 0) {
    stop("bplcutfit: lambda <= 0")
  }
  
  # Log-likelihood positiva
  if (!is.finite(result$loglik)) {
    stop("bplcutfit: log-likelihood not ended")
  }
  
  cat("bplcutfit test ended\n\n")
  return(TRUE)
}


test_bplcutfit()
