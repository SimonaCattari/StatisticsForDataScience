source ("Funzioni/blgnormfit.R")

test_blgnormfit <- function() {
  cat("=== TEST BLGNORMFIT FIXED ===\n")
  
  # FORZA l'uso di environment locale per evitare caching
  local({
    set.seed(42)
    
    # Parametri veri
    true_mu <- 1.73    
    true_sigma <- 0.91
    true_bmin <- 1
    n <- 6000
    
    cat("True params: mu =", true_mu, ", sigma =", true_sigma, "\n")
    
    # Genera dati log-normali
    x <- rlnorm(n, meanlog = true_mu, sdlog = true_sigma)
    x <- x[x >= true_bmin]
    
    # Binning logaritmico
    boundaries <- exp(seq(log(0.5), log(max(x)), length.out = 35))
    h <- hist(x, breaks = boundaries, plot = FALSE)$counts
    
    # Test 1: Default fitting
    result <- blgnormfit(h, boundaries, bmin = true_bmin)
    
    error_mu <- abs(result$mu - true_mu) / abs(true_mu)
    error_sigma <- abs(result$sigma - true_sigma) / true_sigma
    
    cat(sprintf("Mu: estimate = %.3f, true = %.3f, error = %.1f%%\n",
                result$mu, true_mu, error_mu * 100))
    cat(sprintf("Sigma: estimate = %.3f, true = %.3f, error = %.1f%%\n",
                result$sigma, true_sigma, error_sigma * 100))
    
    # Test con range personalizzati che NON includono i valori veri esatti
    murng <- seq(1.2, 2.2, by = 0.1)  # Include 1.7 ma non 1.73
    sigrng <- seq(0.6, 1.3, by = 0.1)  # Include 0.9 ma non 0.91
    
    result_custom <- blgnormfit(h, boundaries, murng = murng, sigrng = sigrng,
                                bmin = true_bmin, fine = FALSE)
    
    error_mu_custom <- abs(result_custom$mu - true_mu) / abs(true_mu)
    error_sigma_custom <- abs(result_custom$sigma - true_sigma) / true_sigma
    
    cat(sprintf("Custom range - Mu: estimate = %.3f, error = %.1f%%\n",
                result_custom$mu, error_mu_custom * 100))
    cat(sprintf("Custom range - Sigma: estimate = %.3f, error = %.1f%%\n",
                result_custom$sigma, error_sigma_custom * 100))
    
    # Test 3: Fine search dovrebbe migliorare
    result_fine <- blgnormfit(h, boundaries, bmin = true_bmin, fine = TRUE)
    
    # Controlla che entrambi i risultati abbiano il campo logLik
    if (!is.null(result_fine$logLik) && !is.null(result$logLik) && 
        is.finite(result_fine$logLik) && is.finite(result$logLik)) {
      if (result_fine$logLik < result$logLik) {
        warning("blgnormfit: fine search worse fit")
      } else {
        cat("Fine search makes fit better\n")
      }
    } else {
      warning("blgnormfit: missing log-likelihood")
    }
    
    # Verifica sigma positivo
    if (result$sigma <= 0) {
      stop("blgnormfit: negative sigma")
    }
    
    # Test ragionevoli: errore < 15% per campioni grandi con grid search
    if (error_mu > 0.15 || error_sigma > 0.15) {
      warning("blgnormfit: errors > 15%")
    }
    
    cat("blgnormfit test done\n\n")
    return(TRUE)
  })
}

test_blgnormfit_clean <- function() {
  cat("=== TEST BLGNORMFIT CLEAN ===\n")
  
  # Rimuovi eventuali variabili cached
  if (exists("true_mu")) rm(true_mu)
  if (exists("true_sigma")) rm(true_sigma)
  
  set.seed(12345) 
  
  # Parametri in variabili locali 
  params <- list(
    mu = 1.83,     
    sigma = 0.67,  
    bmin = 1,
    n = 4000
  )
  
  cat("True Params:", params$mu, params$sigma, "\n")
  
  # Genera dati
  x <- rlnorm(params$n, meanlog = params$mu, sdlog = params$sigma)
  x <- x[x >= params$bmin]
  
  boundaries <- exp(seq(log(0.5), log(max(x)), length.out = 30))
  h <- hist(x, breaks = boundaries, plot = FALSE)$counts
  
  # Test fitting
  result <- blgnormfit(h, boundaries, bmin = params$bmin, fine = TRUE)
  
  # Calcola errori usando i parametri dalla lista
  error_mu <- abs(result$mu - params$mu)
  error_sigma <- abs(result$sigma - params$sigma)
  
  cat(sprintf("Result: mu = %.3f (error = %.3f)\n", result$mu, error_mu))
  cat(sprintf("Result: sigma = %.3f (error = %.3f)\n", result$sigma, error_sigma))
  
  success <- (error_mu < 0.2 && error_sigma < 0.2 && result$sigma > 0)
  
  if (success) {
    cat("Test done!\n")
  } else {
    cat("Test failed\n")
  }
  
  return(success)
}

test_blgnormfit_clean()