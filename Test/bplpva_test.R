# ========================================================================
# TEST
# ========================================================================

source("Funzioni/bplpva.R")

test_bplpva_fixed <- function() {
  cat("=== Test bplpva ===\n\n")
  
  true_alpha <- 2.5
  true_bmin <- 16
  N <- 10000
  
  set.seed(42)
  
  # =====================================================================
  # TEST 1: POWER LAW VERA - P-VALUES DOVREBBERO ESSERE ~UNIFORM(0,1)
  # =====================================================================
  cat("TEST 1: Power Law Vera (p-values alti)\n")
  cat(paste(rep("-", 50), collapse=""), "\n")
  
  p_values <- numeric(20)  # 20 repliche per test statistico
  
  for (i in 1:20) {
    set.seed(10000 + i)
    
    # Genera power law vera
    u <- runif(N)
    x <- true_bmin * (1 - u)^(-1/(true_alpha - 1))
    
    # Binning logaritmico
    boundaries <- true_bmin * 2^(0:10)
    if (max(x) >= max(boundaries)) {
      boundaries[length(boundaries)] <- max(x) * 1.1
    }
    
    h <- hist(x, breaks = boundaries, plot = FALSE, right = FALSE)$counts
    
    # Test 
    result <- bplpva(h, boundaries, true_bmin, true_alpha, 
                     reps = 500, silent = TRUE, seed = i)
    
    p_values[i] <- result$p
    
    if (i <= 5) {
      cat(sprintf("Replica %d: p = %.4f, D* = %.6f\n", i, result$p, result$Dstar))
    }
  }
  
  # Analisi statistiche
  valid_p <- p_values[!is.na(p_values)]
  mean_p <- mean(valid_p)
  acceptance_rate <- sum(valid_p >= 0.1) / length(valid_p)
  
  cat(sprintf("\nRISULTATI TEST 1:\n"))
  cat(sprintf("P-values validi: %d/20\n", length(valid_p)))
  cat(sprintf("Media p-value: %.3f (target ≈ 0.5)\n", mean_p))
  cat(sprintf("Accettazione (p≥0.1): %.0f%% (target ≥ 90%%)\n", 100*acceptance_rate))
  
  # Test distribuzione uniform
  # K-S per testare uniformità
  if (length(valid_p) >= 10) {
    ks_test <- ks.test(valid_p, punif)
    cat(sprintf("Test uniformità: p = %.3f (>0.05 = uniforme)\n", ks_test$p.value))
  }
  
  # =====================================================================
  # TEST 2: LOG-NORMAL - DOVREBBE ESSERE RIGETTATA PER N GRANDI
  # =====================================================================
  cat(sprintf("\nTEST 2: Log-normal (rigetta)\n"))
  cat(paste(rep("-", 50), collapse=""), "\n")
  
  # Parametri log-normal dal paper
  mu_ln <- 0.3
  sigma_ln <- 2.0
  
  p_values_ln <- numeric(10)
  
  for (i in 1:10) {
    set.seed(20000 + i)
    
    # Genera log-normal
    x_ln <- rlnorm(N, meanlog = mu_ln, sdlog = sigma_ln)
    x_ln <- x_ln[x_ln >= true_bmin]  # Tronca sotto bmin
    
    if (length(x_ln) < N*0.8) next  # Skip se troppi dati persi
    
    # Binning
    boundaries_ln <- true_bmin * 2^(0:10)
    if (max(x_ln) >= max(boundaries_ln)) {
      boundaries_ln[length(boundaries_ln)] <- max(x_ln) * 1.1
    }
    
    h_ln <- hist(x_ln, breaks = boundaries_ln, plot = FALSE, right = FALSE)$counts
    
    # Test power law su dati log-normal
    result_ln <- bplpva(h_ln, boundaries_ln, true_bmin, true_alpha, 
                        reps = 500, silent = TRUE, seed = i)
    
    p_values_ln[i] <- result_ln$p
    
    if (i <= 3) {
      cat(sprintf("Replica %d: p = %.4f (basso)\n", i, result_ln$p))
    }
  }
  
  valid_p_ln <- p_values_ln[!is.na(p_values_ln)]
  mean_p_ln <- mean(valid_p_ln)
  rejection_rate <- sum(valid_p_ln < 0.1) / length(valid_p_ln)
  
  cat(sprintf("\nRISULTATI TEST 2:\n"))
  cat(sprintf("P-values validi: %d/10\n", length(valid_p_ln)))
  cat(sprintf("Media p-value: %.3f (target < 0.1)\n", mean_p_ln))
  cat(sprintf("Rigetto (p<0.1): %.0f%% (target ≥ 70%%)\n", 100*rejection_rate))
  
  # =====================================================================
  # VALUTAZIONE FINALE
  # =====================================================================
  cat(sprintf("\n"))
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("VALUTAZIONE FINALE\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Criteri di successo
  power_law_ok <- (mean_p >= 0.4) && (acceptance_rate >= 0.85)
  lognormal_ok <- (mean_p_ln <= 0.2) && (rejection_rate >= 0.6)
  
  cat(sprintf("Power Law Test: %s\n", if(power_law_ok) "PASSED" else "FAILED"))
  cat(sprintf("Media p-value: %.3f %s\n", mean_p, if(mean_p >= 0.4) "(OK)" else "(ERROR)"))
  cat(sprintf("Accettazione: %.0f%% %s\n", 100*acceptance_rate, if(acceptance_rate >= 0.85) "(OK)" else "(ERROR)"))
  
  cat(sprintf("\nLog-normal Test: %s\n", if(lognormal_ok) "PASSED" else "FAILED"))
  cat(sprintf("  Media p-value: %.3f %s\n", mean_p_ln, if(mean_p_ln <= 0.2) "(OK)" else "(ERROR)"))
  cat(sprintf("  Rigetto: %.0f%% %s\n", 100*rejection_rate, if(rejection_rate >= 0.6) "(OK)" else "(ERROR)"))
  
  overall_success <- power_law_ok && lognormal_ok
  
  cat(sprintf("\nRISULTATO COMPLESSIVO: %s\n", 
              if(overall_success) "plpva FUNZIONA CORRETTAMENTE!" else "Servono ulteriori correzioni"))
  
  return(list(
    power_law_ok = power_law_ok,
    lognormal_ok = lognormal_ok,
    overall_success = overall_success,
    mean_p_powerlaw = mean_p,
    mean_p_lognormal = mean_p_ln
  ))
}

# ========================================================================
# ESECUZIONE TEST
# ========================================================================

cat("TEST bplpva\n")
cat(paste(rep("=", 60), collapse=""), "\n")

results <- test_bplpva_fixed()

cat(paste(rep("=", 60), collapse=""), "\n")
cat("Test completato! Verifica i risultati sopra.\n")")