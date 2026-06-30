# ===================================================================
# TEST COMPLETO 
# ===================================================================

source ("Funzioni/bplfit.R")

results_global <- tibble(
  scenario    = character(),
  alpha_true  = double(),
  alpha_hat   = double(),
  bmin_true   = double(),
  bmin_hat    = double(),
)

test_bplfit_complete <- function() {
  cat("Testing power-law parameter estimation from binned data\n\n")
  
  validate_bplfit_result <- function(result, expected_alpha, expected_bmin, test_name, tolerance_alpha = 0.1) {
    # Test 1: Struttura output
    if (is.null(result)) {
      warning(sprintf("%s: bplfit returned NULL", test_name))
      return(list(alpha_error = Inf, alpha_rel_error = Inf, valid = FALSE))
    }
    
    required_fields <- c("alpha", "bmin", "logLik")
    missing_fields <- setdiff(required_fields, names(result))
    if (length(missing_fields) > 0) {
      warning(sprintf("%s: Missing fields: %s", test_name, paste(missing_fields, collapse=", ")))
      return(list(alpha_error = Inf, alpha_rel_error = Inf, valid = FALSE))
    }
    
    # Test 2: Controllo validità alpha PRIMA di calcolare errori
    if (is.na(result$alpha) || !is.finite(result$alpha) || result$alpha <= 1) {
      warning(sprintf("%s: Invalid alpha value: %s", test_name, result$alpha))
      return(list(alpha_error = Inf, alpha_rel_error = Inf, valid = FALSE))
    }
    
    # Test 3: Controllo validità bmin
    if (is.na(result$bmin) || !is.finite(result$bmin) || result$bmin <= 0) {
      warning(sprintf("%s: Invalid bmin value: %s", test_name, result$bmin))
      return(list(alpha_error = Inf, alpha_rel_error = Inf, valid = FALSE))
    }
    
    # Test 4: Controllo validità log-likelihood
    if (is.na(result$logLik) || !is.finite(result$logLik)) {
      warning(sprintf("%s: Invalid log-likelihood: %s", test_name, result$logLik))
    }
    
    # Test 5: Alpha accuracy 
    alpha_error <- abs(result$alpha - expected_alpha)
    alpha_rel_error <- alpha_error / expected_alpha
    
    if (alpha_error > tolerance_alpha) {
      warning(sprintf("%s: Large alpha error. Expected: %.3f, Got: %.3f, Error: %.3f", 
                      test_name, expected_alpha, result$alpha, alpha_error))
    }
    
    # Test 6: bmin plausibility 
    bmin_error <- abs(result$bmin - expected_bmin)
    
    cat(sprintf("%s: α=%.3f (err=%.3f), bmin=%.1f (err=%.1f), logLik=%.2f\n", 
                test_name, result$alpha, alpha_error, result$bmin, bmin_error, result$logLik))
    
    # salva su DF globale
    results_global <<- results_global %>%
      bind_rows(tibble(
        scenario   = test_name,
        alpha_true = expected_alpha,
        alpha_hat  = result$alpha,
        bmin_true  = expected_bmin,
        bmin_hat   = result$bmin,
        bmin_error = abs(result$bmin - expected_bmin)
      ))
    
    list(alpha_error = alpha_error,
         bmin_error  = bmin_error,
         valid       = alpha_error < tolerance_alpha)
  }
  
  # ===================================================================
  # GENERAZIONE DATI POWER-LAW
  # ===================================================================
  generate_powerlaw_data <- function(N, alpha, xmin) {
    if (alpha <= 1) stop("Alpha must be > 1")
    if (N <= 0) stop("N must be positive")
    if (xmin <= 0) stop("xmin must be positive")
    
    u <- runif(N)
    x <- xmin * (1 - u)^(-1 / (alpha - 1))
    return(x)
  }
  
  # ===================================================================
  # BINNING
  # ===================================================================
  safe_binning <- function(x, boundaries, test_name) {
    # Estendi boundaries se necessario
    if (min(x) < min(boundaries)) {
      boundaries <- c(min(x) * 0.9, boundaries)
    }
    if (max(x) >= max(boundaries)) {
      boundaries[length(boundaries)] <- max(x) * 1.1
    }
    
    h <- tryCatch({
      hist(x, breaks = boundaries, plot = FALSE, right = FALSE)$counts
    }, error = function(e) {
      warning(sprintf("%s: Binning failed: %s", test_name, e$message))
      return(NULL)
    })
    
    return(list(h = h, boundaries = boundaries))
  }
  
  # ===================================================================
  # TEST 1
  # ===================================================================
  cat("Test 1: Base Case - Standard Parameters\n")
  
  set.seed(42)  
  
  N <- 10000
  xmin_true <- 10
  alpha_true <- 2.5
  
  # Generazione dati power-law
  x <- generate_powerlaw_data(N, alpha_true, xmin_true)
  
  # Schema binning logaritmico
  boundaries <- 10 * 2^(0:15)
  
  # Binning 
  binning_result <- safe_binning(x, boundaries, "Base-Case")
  if (is.null(binning_result$h)) {
    cat("Base-Case: Binning failed\n")
  } else {
    # Chiamata a bplfit
    result_base <- tryCatch({
      bplfit(binning_result$h, binning_result$boundaries)
    }, error = function(e) {
      warning(sprintf("Base-Case: bplfit failed: %s", e$message))
      return(NULL)
    })
    
    errors_base <- validate_bplfit_result(result_base, alpha_true, xmin_true, "Base-Case")
  }
  
  # ===================================================================
  # TEST 2: DIVERSI VALORI DI ALPHA
  # ===================================================================
  cat("\nTest 2: Different Alpha Values\n")
  
  alpha_values <- c(1.5, 2.0, 2.5, 3.0, 3.5)
  alpha_results <- list()
  
  for (i in seq_along(alpha_values)) {
    set.seed(1000 + i)
    alpha_test <- alpha_values[i]
    
    # Generazione dati
    x_test <- generate_powerlaw_data(N, alpha_test, xmin_true)
    
    # Binning
    binning_test <- safe_binning(x_test, boundaries, sprintf("Alpha-%.1f", alpha_test))
    
    if (!is.null(binning_test$h)) {
      # Fit
      result_test <- tryCatch({
        bplfit(binning_test$h, binning_test$boundaries)
      }, error = function(e) {
        warning(sprintf("Alpha-%.1f: bplfit failed: %s", alpha_test, e$message))
        return(NULL)
      })
      
      errors_test <- validate_bplfit_result(result_test, alpha_test, xmin_true, 
                                            sprintf("Alpha-%.1f", alpha_test))
      alpha_results[[i]] <- errors_test
    } else {
      alpha_results[[i]] <- list(alpha_error = Inf, valid = FALSE)
    }
  }
  
  # Statistiche sugli errori
  valid_errors <- sapply(alpha_results, function(x) if(x$valid) x$alpha_error else NA)
  valid_errors <- valid_errors[!is.na(valid_errors) & is.finite(valid_errors)]
  
  if (length(valid_errors) > 0) {
    mean_alpha_error <- mean(valid_errors)
    max_alpha_error <- max(valid_errors)
    cat(sprintf("  Alpha errors: mean=%.4f, max=%.4f (from %d valid tests)\n", 
                mean_alpha_error, max_alpha_error, length(valid_errors)))
    
    if (mean_alpha_error > 0.1) {
      warning("High mean alpha error across different alpha values")
    }
  } else {
    cat("No valid alpha estimates obtained\n")
    mean_alpha_error <- Inf
  }
  
  # ===================================================================
  # TEST 3: DIVERSI SCHEMI DI BINNING
  # ===================================================================
  cat("\nTest 3: Different Binning Schemes\n")
  
  set.seed(42)
  x_binning <- generate_powerlaw_data(N, alpha_true, xmin_true)
  
  # Test A: Logaritmico fine (c=2)
  boundaries_fine <- 10 * 2^(seq(0, log2(max(x_binning)/10) + 1, by = 0.5))
  binning_fine <- safe_binning(x_binning, boundaries_fine, "Fine-Binning")
  
  if (!is.null(binning_fine$h)) {
    result_fine <- tryCatch({
      bplfit(binning_fine$h, binning_fine$boundaries)
    }, error = function(e) NULL)
    validate_bplfit_result(result_fine, alpha_true, xmin_true, "Fine-Binning")
  }
  
  # Test B: Logaritmico grosso (c=4)  
  boundaries_coarse <- 10 * 4^(0:8)
  binning_coarse <- safe_binning(x_binning, boundaries_coarse, "Coarse-Binning")
  
  if (!is.null(binning_coarse$h)) {
    result_coarse <- tryCatch({
      bplfit(binning_coarse$h, binning_coarse$boundaries)
    }, error = function(e) NULL)
    validate_bplfit_result(result_coarse, alpha_true, xmin_true, "Coarse-Binning", tolerance_alpha = 0.15)
  }
  
  # Test C: Lineare (per confronto)
  boundaries_linear <- seq(10, ceiling(max(x_binning)), by = 50)
  binning_linear <- safe_binning(x_binning, boundaries_linear, "Linear-Binning")
  
  if (!is.null(binning_linear$h)) {
    result_linear <- tryCatch({
      bplfit(binning_linear$h, binning_linear$boundaries)
    }, error = function(e) NULL)
    validate_bplfit_result(result_linear, alpha_true, xmin_true, "Linear-Binning", tolerance_alpha = 0.3)
  }
  
  # Confronto qualità fit
  if (exists("result_fine") && exists("result_coarse") && exists("result_linear") &&
      !is.null(result_fine) && !is.null(result_coarse) && !is.null(result_linear)) {
    cat(sprintf("Binning comparison - LogLik: Fine=%.2f, Coarse=%.2f, Linear=%.2f\n",
                result_fine$logLik, result_coarse$logLik, result_linear$logLik))
  }
  
  # ===================================================================
  # TEST 4: EFFETTO DIMENSIONE CAMPIONE
  # ===================================================================
  cat("\nTest 4: Sample Size Effect\n")
  
  sample_sizes <- c(1000, 5000, 10000, 50000)
  size_results <- list()
  
  for (i in seq_along(sample_sizes)) {
    set.seed(3000 + i)
    N_test <- sample_sizes[i]
    
    # Generazione dati
    x_size <- generate_powerlaw_data(N_test, alpha_true, xmin_true)
    
    # Binning
    binning_size <- safe_binning(x_size, boundaries, sprintf("N-%d", N_test))
    
    if (!is.null(binning_size$h)) {
      # Fit
      result_size <- tryCatch({
        bplfit(binning_size$h, binning_size$boundaries)
      }, error = function(e) NULL)
      
      errors_size <- validate_bplfit_result(result_size, alpha_true, xmin_true, 
                                            sprintf("N-%d", N_test))
      size_results[[i]] <- errors_size
    } else {
      size_results[[i]] <- list(alpha_error = Inf, valid = FALSE)
    }
  }
  
  valid_size_errors <- sapply(size_results, function(x) if(x$valid) x$alpha_error else NA)
  valid_indices <- which(!is.na(valid_size_errors) & is.finite(valid_size_errors))
  
  if (length(valid_indices) >= 2) {
    first_error <- valid_size_errors[valid_indices[1]]
    last_error <- valid_size_errors[valid_indices[length(valid_indices)]]
    cat(sprintf("Sample size trend: N=%dk err=%.4f, N=%dk err=%.4f\n", 
                sample_sizes[valid_indices[1]]/1000, first_error,
                sample_sizes[valid_indices[length(valid_indices)]]/1000, last_error))
    
    if (last_error > first_error * 1.5) {
      warning("Alpha error should decrease with larger sample size")
    }
  }
  
  # ===================================================================
  # TEST 5: DATI CON RUMORE
  # ===================================================================
  cat("\nTest 5: Robustness with Noisy Data\n")
  
  set.seed(42)
  
  # Genera distribuzione composita (power-law + noise)
  N_pl <- 8000  # 80% power-law
  N_noise <- 2000  # 20% noise uniforme
  
  x_pl <- generate_powerlaw_data(N_pl, alpha_true, xmin_true)
  x_noise <- runif(N_noise, min = 1, max = xmin_true)  # Noise sotto bmin
  
  x_composite <- c(x_pl, x_noise)
  
  # Binning
  binning_robust <- safe_binning(x_composite, boundaries, "Noisy-Data")
  
  if (!is.null(binning_robust$h)) {
    # Fit (dovrebbe trovare la regione power-law automaticamente)
    result_robust <- tryCatch({
      bplfit(binning_robust$h, binning_robust$boundaries)
    }, error = function(e) NULL)
    
    errors_robust <- validate_bplfit_result(result_robust, alpha_true, xmin_true, 
                                            "Noisy-Data", tolerance_alpha = 0.2)
    
    # Test: bmin esclude il noise
    if (!is.null(result_robust) && errors_robust$valid && result_robust$bmin < xmin_true) {
      warning("Robust test: bmin should be >= xmin_true to exclude noise")
    }
  }
  
  # ===================================================================
  # TEST 6: CONFRONTO CON METODI NAIVI
  # ===================================================================
  cat("\nTest 6: Comparison with Naive Methods\n")
  
  set.seed(42)
  x_compare <- generate_powerlaw_data(N, alpha_true, xmin_true)
  
  # Binning per confronto
  binning_compare <- safe_binning(x_compare, boundaries, "MLE-vs-Naive")
  
  if (!is.null(binning_compare$h)) {
    # Metodo MLE (bplfit)
    result_mle <- tryCatch({
      bplfit(binning_compare$h, binning_compare$boundaries)
    }, error = function(e) NULL)
    
    if (!is.null(result_mle) && !is.na(result_mle$alpha)) {
      # Metodo naive: regressione lineare su log-log
      bin_centers <- sqrt(binning_compare$boundaries[-length(binning_compare$boundaries)] * 
                            binning_compare$boundaries[-1])
      valid_bins <- binning_compare$h > 0
      
      if (sum(valid_bins) >= 3) {
        log_centers <- log(bin_centers[valid_bins])
        log_freqs <- log(binning_compare$h[valid_bins])
        
        # Regressione lineare
        lm_fit <- tryCatch({
          lm(log_freqs ~ log_centers)
        }, error = function(e) NULL)
        
        if (!is.null(lm_fit)) {
          alpha_naive <- -coef(lm_fit)[2]
          
          # Confronto (solo se entrambi validi)
          if (is.finite(alpha_naive) && alpha_naive > 1) {
            alpha_mle_error <- abs(result_mle$alpha - alpha_true)
            alpha_naive_error <- abs(alpha_naive - alpha_true)
            
            cat(sprintf("Method comparison:\n"))
            cat(sprintf("MLE: α=%.3f (error=%.4f)\n", result_mle$alpha, alpha_mle_error))
            cat(sprintf("Naive: α=%.3f (error=%.4f)\n", alpha_naive, alpha_naive_error))
            
            if (alpha_naive_error < alpha_mle_error) {
              warning("Unexpected: naive method more accurate than MLE")
            }
          }
        }
      }
    }
  }
  
  # ===================================================================
  # TEST 7: EDGE CASES SEMPLIFICATI
  # ===================================================================
  cat("\nTest 7: Edge Cases\n")
  
  # Test A: Pochi bin nella coda
  set.seed(42)
  x_few <- generate_powerlaw_data(1000, alpha_true, xmin_true)
  
  boundaries_few <- c(5, 10, 20, 50, 100, max(x_few) * 1.1)
  binning_few <- safe_binning(x_few, boundaries_few, "Few-Bins")
  
  if (!is.null(binning_few$h)) {
    result_few <- tryCatch({
      bplfit(binning_few$h, binning_few$boundaries)
    }, error = function(e) {
      cat(sprintf("Few-bins case handled gracefully: %s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(result_few)) {
      validate_bplfit_result(result_few, alpha_true, xmin_true, "Few-Bins", tolerance_alpha = 0.5)
    }
  }
  
  # Test B: Alpha estremi
  set.seed(7000)
  alpha_extreme <- 1.1  # Molto vicino al limite teorico
  
  x_extreme <- tryCatch({
    generate_powerlaw_data(N, alpha_extreme, xmin_true)
  }, error = function(e) {
    cat(sprintf("Extreme alpha correctly rejected: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(x_extreme)) {
    binning_extreme <- safe_binning(x_extreme, boundaries, "Extreme-Alpha")
    
    if (!is.null(binning_extreme$h)) {
      result_extreme <- tryCatch({
        bplfit(binning_extreme$h, binning_extreme$boundaries)
      }, error = function(e) {
        cat(sprintf("Extreme alpha handled: %s\n", e$message))
        return(NULL)
      })
      
      if (!is.null(result_extreme)) {
        validate_bplfit_result(result_extreme, alpha_extreme, xmin_true, "Extreme-Alpha", tolerance_alpha = 0.3)
      }
    }
  }
  
  return(list(
    base_result = if(exists("result_base")) result_base else NULL,
    mean_alpha_error = if(exists("mean_alpha_error")) mean_alpha_error else Inf,
    tests_completed = TRUE
  ))
}

# ===================================================================
# FUNZIONE DI TEST RAPIDO CORRETTA
# ===================================================================

test_bplfit_quick <- function() {
  cat("=== QUICK TEST bplfit ===\n")

  set.seed(42)

  # Parametri base
  N <- 5000
  xmin_true <- 10
  alpha_true <- 2.5

  # Genera e binna
  x <- xmin_true * (runif(N))^(-1 / (alpha_true - 1))
  boundaries <- 10 * 2^(0:12)
  if (max(x) >= max(boundaries)) {
    boundaries[length(boundaries)] <- max(x) * 1.1
  }
  h <- hist(x, breaks = boundaries, plot = FALSE, right = FALSE)$counts

  # Test
  result <- tryCatch({
    bplfit(h, boundaries)
  }, error = function(e) {
    cat("Quick test FAILED: bplfit error\n")
    cat(sprintf("Error: %s\n", e$message))
    return(NULL)
  })

  # CONTROLLO VALIDITÀ RISULTATO
  if (is.null(result)) {
    cat("Quick test FAILED: bplfit returned NULL\n")
    return(FALSE)
  }

  if (is.na(result$alpha) || is.na(result$bmin) || !is.finite(result$alpha) || !is.finite(result$bmin)) {
    cat("Quick test FAILED: Invalid result from bplfit\n")
    cat(sprintf("  result: alpha=%s, bmin=%s, logLik=%s\n",
                result$alpha, result$bmin, result$logLik))
    return(FALSE)
  }

  # Validazione veloce
  alpha_error <- abs(result$alpha - alpha_true)
  bmin_error <- abs(result$bmin - xmin_true)

  alpha_ok <- alpha_error < 0.15
  bmin_ok <- bmin_error < 20  # Tolleranza larga per bmin
  loglik_ok <- is.finite(result$logLik)

  if (alpha_ok && bmin_ok && loglik_ok) {
    cat("Quick test PASSED\n")
    cat(sprintf("α=%.3f (target=%.1f, err=%.3f), bmin=%.1f (target=%d, err=%.1f)\n",
                result$alpha, alpha_true, alpha_error, result$bmin, xmin_true, bmin_error))
    return(TRUE)
  } else {
    cat("Quick test FAILED\n")
    cat(sprintf("α=%.3f (err=%.3f, ok=%s), bmin=%.1f (err=%.1f, ok=%s), logLik finite=%s\n",
                result$alpha, alpha_error, alpha_ok,
                result$bmin, bmin_error, bmin_ok, loglik_ok))
    return(FALSE)
  }
}

# ===================================================================
# TEST DI DEBUG MANUALE
# ===================================================================

test_bplfit_debug <- function() {
  cat("=== DEBUG TEST bplfit ===\n")

  # Test manuale semplificato per isolare problemi
  boundaries <- c(1, 10, 100)  # Solo 2 bin
  h <- c(90, 10)               # 90% primo bin, 10% secondo

  cat("Debug input:\n")
  cat(sprintf("boundaries: %s\n", paste(boundaries, collapse=", ")))
  cat(sprintf("h: %s\n", paste(h, collapse=", ")))

  result <- tryCatch({
    bplfit(h, boundaries)
  }, error = function(e) {
    cat(sprintf("Debug FAILED: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(result)) {
    cat("Debug result:\n")
    cat(sprintf("alpha: %s\n", result$alpha))
    cat(sprintf("bmin: %s\n", result$bmin))
    cat(sprintf("logLik: %s\n", result$logLik))
    cat(sprintf("D: %s\n", result$D))
  }

  return(result)
}

# ===================================================================
# ESECUZIONE PRINCIPALE
# ===================================================================

if (interactive()) {
  cat("Running comprehensive bplfit tests...\n\n")

  # Test rapido prima
  quick_success <- test_bplfit_quick() 

  cat("\n", paste(rep("=", 50), collapse=""), "\n")

  if (quick_success) {
    # Test completo solo se quick test passa
    results <- test_bplfit_complete()

    cat("\nAll bplfit tests completed!\n")
    if (is.finite(results$mean_alpha_error)) {
      cat(sprintf("Overall alpha accuracy: %.4f average error\n", results$mean_alpha_error))
    } else {
      cat("Some tests failed - check warnings above\n")
    }
  } else {
    cat("\n Failed\n")
    test_bplfit_quick()
  }

}