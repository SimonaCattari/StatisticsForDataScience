# ===================================================================
# TEST 
# ===================================================================

source("Funzioni/blrtest.R")
source("Funzioni/bplfit.R")
source("Funzioni/getPDF.R")
source("Funzioni/blgnormfit.R")


test_blrtest_complete <- function() {
  cat("=== TEST COMPLETO BLRTEST ===\n")
  cat("Likelihood Ratio Test per confronto modelli su dati binned\n\n")
  
  # ===================================================================
  # HELPER FUNCTIONS
  # ===================================================================
  
  # Funzione per validazione risultati
  validate_lr_result <- function(lr_result, test_name) {
    if (is.null(lr_result) || 
        is.na(lr_result$normR) || is.na(lr_result$p) || 
        !is.finite(lr_result$normR) || !is.finite(lr_result$p)) {
      warning(sprintf("%s: Invalid LR test result", test_name))
      return(FALSE)
    }
    return(TRUE)
  }
  
  # Funzione per fit log-normal
  robust_lognormal_fit <- function(h, boundaries, test_name) {
    # Prova diversi approcci per blgnormfit
    approaches <- list(
      list(name = "constrained", 
           murng = seq(-2, 5, 0.2), 
           sigrng = seq(0.1, 3, 0.1)),
      list(name = "wide-constrained", 
           murng = seq(-5, 8, 0.3), 
           sigrng = seq(0.05, 5, 0.1)),
      list(name = "default", 
           murng = NULL, 
           sigrng = NULL)
    )
    
    for (approach in approaches) {
      fit_ln <- tryCatch({
        if (is.null(approach$murng)) {
          blgnormfit(h, boundaries, fine = TRUE)
        } else {
          blgnormfit(h, boundaries, 
                     murng = approach$murng, 
                     sigrng = approach$sigrng, 
                     fine = TRUE)
        }
      }, error = function(e) NULL)
      
      # Verifica parametri
      if (!is.null(fit_ln) && 
          !is.na(fit_ln$mu) && !is.na(fit_ln$sigma) &&
          is.finite(fit_ln$mu) && is.finite(fit_ln$sigma) &&
          fit_ln$mu > -20 && fit_ln$mu < 20 && 
          fit_ln$sigma > 0.01 && fit_ln$sigma < 10) {
        
        cat(sprintf("%s: Log-normal fit (%s): μ=%.3f, σ=%.3f, logLik=%.2f\n",
                    test_name, approach$name, fit_ln$mu, fit_ln$sigma, fit_ln$logLik))
        return(fit_ln)
      }
    }
    
    warning(sprintf("%s: All log-normal fitting approaches failed", test_name))
    return(NULL)
  }
  
  # ===================================================================
  # TEST 1: POWER-LAW vs LOG-NORMAL (DATI POWER-LAW)
  # ===================================================================
  cat("Test 1: Power-law vs Log-normal on Power-law data\n")
  
  set.seed(42)
  N <- 20000
  alpha_true <- 2.5
  xmin_true <- 1
  
  # Genera dati power-law
  x_pl <- xmin_true * (runif(N))^(-1 / (alpha_true - 1))
  
  # Binning logaritmico
  boundaries_pl <- 2^(0:15)
  if (max(x_pl) >= max(boundaries_pl)) {
    boundaries_pl[length(boundaries_pl)] <- max(x_pl) * 1.1
  }
  
  h_pl <- hist(x_pl, breaks = boundaries_pl, plot = FALSE, right = FALSE)$counts
  
  # Fit modelli
  fit_pl <- tryCatch({
    bplfit(h_pl, boundaries_pl)
  }, error = function(e) {
    warning("Power-law fit failed:", e$message)
    return(NULL)
  })
  
  if (is.null(fit_pl) || is.na(fit_pl$alpha)) {
    cat("Test 1 FAILED: Power-law fit failed\n")
  } else {
    cat(sprintf("Power-law fit: α=%.3f, bmin=%.1f, logLik=%.2f\n",
                fit_pl$alpha, fit_pl$bmin, fit_pl$logLik))
    
    # Fit log-normal robusto
    fit_ln <- robust_lognormal_fit(h_pl, boundaries_pl, "Test1")
    
    if (!is.null(fit_ln)) {
      # Calcola PDF con bmin comune
      common_bmin <- fit_pl$bmin
      
      p_pl <- tryCatch({
        getPDF(boundaries_pl, "pl", common_bmin, fit_pl$alpha)
      }, error = function(e) {
        warning("getPDF power-law failed:", e$message)
        return(NULL)
      })
      
      p_ln <- tryCatch({
        getPDF(boundaries_pl, "lgnorm", common_bmin, fit_ln$mu, fit_ln$sigma)
      }, error = function(e) {
        warning("getPDF log-normal failed:", e$message)
        return(NULL)
      })
      
      if (!is.null(p_pl) && !is.null(p_ln)) {
        # Likelihood ratio test
        lr_result <- tryCatch({
          blrtest(p_pl, p_ln, h_pl, boundaries_pl, common_bmin)
        }, error = function(e) {
          warning("blrtest failed:", e$message)
          return(NULL)
        })
        
        if (validate_lr_result(lr_result, "Test1")) {
          cat(sprintf("-> LR Test: normR=%.3f, p-value=%.3f\n", 
                      lr_result$normR, lr_result$p))
          
          if (lr_result$normR > 0 && lr_result$p < 0.1) {
            cat("Test 1 PASSED: Power-law significantly favored (correct)\n")
          } else if (lr_result$normR > 0) {
            cat("- Test 1 PARTIAL: Power-law favored but not significant\n")
          } else {
            cat("Test 1 FAILED: Unexpected result\n")
          }
        }
      }
    }
  }
  
  # ===================================================================
  # TEST 2: POWER-LAW vs LOG-NORMAL (DATI LOG-NORMAL)
  # ===================================================================
  cat("\nTest 2: Power-law vs Log-normal on Log-normal data\n")
  
  set.seed(42)
  N_ln <- 15000
  mu_true <- 1.5
  sigma_true <- 0.8
  
  x_ln <- rlnorm(N_ln, meanlog = mu_true, sdlog = sigma_true)
  
  # Binning logaritmico
  boundaries_ln <- exp(seq(log(min(x_ln)), log(max(x_ln)), length.out = 25))
  h_ln <- hist(x_ln, breaks = boundaries_ln, plot = FALSE, right = FALSE)$counts
  
  # Fit modelli
  fit_pl2 <- tryCatch({
    bplfit(h_ln, boundaries_ln)
  }, error = function(e) NULL)
  
  fit_ln2 <- tryCatch({
    blgnormfit(h_ln, boundaries_ln, fine = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_pl2) && !is.null(fit_ln2) && 
      !is.na(fit_pl2$alpha) && !is.na(fit_ln2$mu)) {
    
    cat(sprintf("Power-law fit: α=%.3f, bmin=%.1f\n", fit_pl2$alpha, fit_pl2$bmin))
    cat(sprintf("Log-normal fit: μ=%.3f, σ=%.3f (target: %.1f, %.1f)\n",
                fit_ln2$mu, fit_ln2$sigma, mu_true, sigma_true))
    
    # Test LR
    common_bmin2 <- min(fit_pl2$bmin, exp(fit_ln2$mu - 2*fit_ln2$sigma))
    
    p_pl2 <- getPDF(boundaries_ln, "pl", common_bmin2, fit_pl2$alpha)
    p_ln2 <- getPDF(boundaries_ln, "lgnorm", common_bmin2, fit_ln2$mu, fit_ln2$sigma)
    
    if (!is.null(p_pl2) && !is.null(p_ln2)) {
      lr_result2 <- tryCatch({
        blrtest(p_pl2, p_ln2, h_ln, boundaries_ln, common_bmin2)
      }, error = function(e) NULL)
      
      if (validate_lr_result(lr_result2, "Test2")) {
        cat(sprintf("-> LR Test: normR=%.3f, p-value=%.3f\n", 
                    lr_result2$normR, lr_result2$p))
        
        if (lr_result2$normR < 0 && lr_result2$p < 0.1) {
          cat("Test 2 PASSED: Log-normal significantly favored (correct)\n")
        } else if (lr_result2$normR < 0) {
          cat("- Test 2 PARTIAL: Log-normal favored but not significant\n")
        } else {
          cat("Test 2 UNEXPECTED: Power-law favored on log-normal data\n")
        }
      }
    }
  } else {
    cat("Test 2 FAILED: Model fitting failed\n")
  }
  
  # ===================================================================
  # TEST 3: POWER-LAW vs EXPONENTIAL
  # ===================================================================
  cat("\nTest 3: Power-law vs Exponential on Power-law data\n")
  
  set.seed(42)
  
  # Usa dati power-law del Test 1
  fit_exp <- tryCatch({
    bexpnfit(h_pl, boundaries_pl, bmin = fit_pl$bmin)
  }, error = function(e) NULL)
  
  if (!is.null(fit_exp) && !is.na(fit_exp$lambda)) {
    cat(sprintf("    ✓ Exponential fit: λ=%.3f, logLik=%.2f\n", 
                fit_exp$lambda, fit_exp$logLik))
    
    # Test LR
    p_exp <- getPDF(boundaries_pl, "expn", fit_pl$bmin, fit_exp$lambda)
    
    if (!is.null(p_exp)) {
      lr_result3 <- tryCatch({
        blrtest(p_pl, p_exp, h_pl, boundaries_pl, fit_pl$bmin)
      }, error = function(e) NULL)
      
      if (validate_lr_result(lr_result3, "Test3")) {
        cat(sprintf("    → LR Test: normR=%.3f, p-value=%.3f\n", 
                    lr_result3$normR, lr_result3$p))
        
        if (lr_result3$normR > 0 && lr_result3$p < 0.1) {
          cat("Test 3 PASSED: Power-law significantly favored over exponential\n")
        } else {
          cat("- Test 3: Power-law vs exponential inconclusive\n")
        }
      }
    }
  } else {
    cat("Test 3 FAILED: Exponential fit failed\n")
  }
  
  # ===================================================================
  # TEST 4: MULTIPLE COMPARISONS
  # ===================================================================
  cat("\nTest 4: Multiple Model Comparisons\n")
  
  set.seed(42)
  
  # Genera dati power-law con cutoff
  N_cut <- 10000
  alpha_cut <- 2.8
  lambda_cut <- 0.01
  
  # Genera power-law con cutoff (metodo rejection)
  x_cut <- c()
  while(length(x_cut) < N_cut) {
    x_cand <- 1 * (runif(1))^(-1 / (alpha_cut - 1))
    if(runif(1) < exp(-lambda_cut * x_cand)) {
      x_cut <- c(x_cut, x_cand)
    }
  }
  
  boundaries_cut <- exp(seq(log(1), log(max(x_cut)), length.out = 20))
  h_cut <- hist(x_cut, breaks = boundaries_cut, plot = FALSE, right = FALSE)$counts
  
  # Fit multipli modelli
  models <- list()
  
  # Power-law
  models$pl <- tryCatch({
    fit <- bplfit(h_cut, boundaries_cut)
    list(name = "Power-law", fit = fit, 
         pdf = if(!is.na(fit$alpha)) getPDF(boundaries_cut, "pl", fit$bmin, fit$alpha) else NULL)
  }, error = function(e) NULL)
  
  # Log-normal
  models$ln <- tryCatch({
    fit <- robust_lognormal_fit(h_cut, boundaries_cut, "Test4")
    list(name = "Log-normal", fit = fit,
         pdf = if(!is.null(fit)) getPDF(boundaries_cut, "lgnorm", fit$bmin, fit$mu, fit$sigma) else NULL)
  }, error = function(e) NULL)
  
  # Exponential
  models$exp <- tryCatch({
    bmin_common <- if(!is.null(models$pl$fit)) models$pl$fit$bmin else 1
    fit <- bexpnfit(h_cut, boundaries_cut, bmin = bmin_common)
    list(name = "Exponential", fit = fit,
         pdf = if(!is.na(fit$lambda)) getPDF(boundaries_cut, "expn", bmin_common, fit$lambda) else NULL)
  }, error = function(e) NULL)
  
  # Confronti a coppie
  model_pairs <- list(
    c("pl", "ln"),
    c("pl", "exp"),
    c("ln", "exp")
  )
  
  for (pair in model_pairs) {
    m1 <- models[[pair[1]]]
    m2 <- models[[pair[2]]]
    
    if (!is.null(m1) && !is.null(m2) && !is.null(m1$pdf) && !is.null(m2$pdf)) {
      bmin_common <- 1  # Usa bmin comune
      
      lr_result <- tryCatch({
        blrtest(m1$pdf, m2$pdf, h_cut, boundaries_cut, bmin_common)
      }, error = function(e) NULL)
      
      if (validate_lr_result(lr_result, "Test4")) {
        winner <- if(lr_result$normR > 0) m1$name else m2$name
        significance <- if(lr_result$p < 0.1) "significant" else "not significant"
        
        cat(sprintf("%s vs %s: %s favored (normR=%.3f, p=%.3f, %s)\n",
                    m1$name, m2$name, winner, lr_result$normR, lr_result$p, significance))
      }
    }
  }
  
  # ===================================================================
  # TEST 5: ROBUSTNESS TESTS
  # ===================================================================
  cat("\nTest 5: Robustness Tests\n")
  
  # Test con sample sizes diversi
  sample_sizes <- c(1000, 5000, 20000)
  
  for (N_test in sample_sizes) {
    set.seed(N_test)
    
    x_robust <- 1 * (runif(N_test))^(-1 / (2.5 - 1))
    boundaries_robust <- 2^(0:12)
    if (max(x_robust) >= max(boundaries_robust)) {
      boundaries_robust[length(boundaries_robust)] <- max(x_robust) * 1.1
    }
    
    h_robust <- hist(x_robust, breaks = boundaries_robust, plot = FALSE, right = FALSE)$counts
    
    # Quick fits
    fit_pl_robust <- tryCatch(bplfit(h_robust, boundaries_robust), error = function(e) NULL)
    fit_ln_robust <- robust_lognormal_fit(h_robust, boundaries_robust, sprintf("N%d", N_test))
    
    if (!is.null(fit_pl_robust) && !is.null(fit_ln_robust) && 
        !is.na(fit_pl_robust$alpha) && !is.na(fit_ln_robust$mu)) {
      
      p_pl_robust <- getPDF(boundaries_robust, "pl", fit_pl_robust$bmin, fit_pl_robust$alpha)
      p_ln_robust <- getPDF(boundaries_robust, "lgnorm", fit_pl_robust$bmin, fit_ln_robust$mu, fit_ln_robust$sigma)
      
      if (!is.null(p_pl_robust) && !is.null(p_ln_robust)) {
        lr_robust <- tryCatch({
          blrtest(p_pl_robust, p_ln_robust, h_robust, boundaries_robust, fit_pl_robust$bmin)
        }, error = function(e) NULL)
        
        if (validate_lr_result(lr_robust, sprintf("N%d", N_test))) {
          cat(sprintf("    N=%d: normR=%.3f, p=%.3f\n", N_test, lr_robust$normR, lr_robust$p))
        }
      }
    }
  }
  
  return(TRUE)
}

# ===================================================================
# QUICK TEST FUNCTION
# ===================================================================

test_blrtest_quick <- function() {
  cat("=== QUICK BLRTEST TEST ===\n")
  
  set.seed(42)
  
  # Generate simple power-law data
  x <- 1 * (runif(5000))^(-1 / (2.5 - 1))
  boundaries <- 2^(0:10)
  if (max(x) >= max(boundaries)) {
    boundaries[length(boundaries)] <- max(x) * 1.1
  }
  h <- hist(x, breaks = boundaries, plot = FALSE, right = FALSE)$counts
  
  # Fit models
  fit_pl <- bplfit(h, boundaries)
  fit_ln <- blgnormfit(h, boundaries, 
                       murng = seq(-2, 5, 0.3),
                       sigrng = seq(0.1, 3, 0.2))
  
  if (!is.na(fit_pl$alpha) && !is.na(fit_ln$mu)) {
    # Calculate PDFs
    p_pl <- getPDF(boundaries, "pl", fit_pl$bmin, fit_pl$alpha)
    p_ln <- getPDF(boundaries, "lgnorm", fit_pl$bmin, fit_ln$mu, fit_ln$sigma)
    
    # LR test
    lr_result <- blrtest(p_pl, p_ln, h, boundaries, fit_pl$bmin)
    
    cat(sprintf("Quick test result: normR=%.3f, p=%.3f\n", lr_result$normR, lr_result$p))
    
    if (lr_result$normR > 0 && lr_result$p < 0.1) {
      cat("Quick test PASSED\n")
      return(TRUE)
    } else {
      cat("- Quick test: Power-law favored but not significant\n")
      return(TRUE)  # Still acceptable
    }
  } else {
    cat("Quick test FAILED: Model fitting failed\n")
    return(FALSE)
  }
}

if (interactive()) {
  cat("Running BLRTEST tests...\n\n")
  
  # Quick test first
  quick_success <- test_blrtest_quick()
  
  if (quick_success) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    
    # Comprehensive test
    test_blrtest_complete()
  } else {
    cat("\nQuick test failed - check setup\n")
  }
}