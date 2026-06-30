# ===================================================================
# TEST COMPLETO
# ===================================================================

source("Funzioni/getPDF.R")  

test_getPDF_complete <- function() {
  cat("=== TEST COMPLETO getPDF ===\n")
  cat("Testing PDF calculation for all distribution types\n\n")
  
  # Helper function per validazione generale
  validate_pdf_output <- function(prden, boundaries, test_name) {
    expected_length <- length(boundaries) - 1
    
    # Test 1: Lunghezza corretta
    if (length(prden) != expected_length) {
      stop(paste(test_name, ": Wrong output length. Expected:", expected_length, "Got:", length(prden)))
    }
    
    # Test 2: Valori non negativi
    if (any(prden < 0)) {
      stop(paste(test_name, ": Negative probabilities found"))
    }
    
    # Test 3: Valori finiti
    if (any(!is.finite(prden))) {
      stop(paste(test_name, ": Non-finite values found"))
    }
    
    # Test 4: Normalizzazione (somma = 1)
    sum_prden <- sum(prden)
    if (abs(sum_prden - 1.0) > 1e-10) {
      stop(paste(test_name, ": Poor normalization. Sum =", sum_prden, "Expected: 1.0"))
    }
    
    cat(sprintf("%s: length=%d, sum=%.10f, range=[%.6f, %.6f]\n", 
                test_name, length(prden), sum_prden, min(prden), max(prden)))
  }
  
  # ===================================================================
  # TEST 1: POWER-LAW
  # ===================================================================
  cat("Test 1: Power-law Distribution\n")
  
  boundaries_pl <- c(1, 10, 100, 1000, Inf)
  bmin_pl <- 1
  alpha_pl <- 2.0
  
  prden_pl <- getPDF(boundaries_pl, "pl", bmin_pl, alpha_pl)
  validate_pdf_output(prden_pl, boundaries_pl, "Power-law")
  
  # Test formula specifica per power-law con normalizzazione
  # P(bi <= X < bi+1) = C/(α-1) * [bi^(1-α) - bi+1^(1-α)]
  l <- boundaries_pl[-length(boundaries_pl)]
  u <- boundaries_pl[-1]
  u[is.infinite(u)] <- 1e10  # Approssimazione per test
  
  # Formula teorica (non normalizzata)
  theoretical_unnorm <- l^(1-alpha_pl) - u^(1-alpha_pl)
  theoretical_norm <- theoretical_unnorm / sum(theoretical_unnorm)
  
  # Confronto (tolleranza per Inf boundary)
  rel_errors <- abs(prden_pl - theoretical_norm) / theoretical_norm
  if (any(rel_errors > 0.01, na.rm = TRUE)) {
    warning("Power-law: Large relative errors in formula verification")
  }
  
  # ===================================================================
  # TEST 2: EXPONENTIAL
  # ===================================================================
  cat("Test 2: Exponential Distribution\n")
  
  boundaries_exp <- c(1, 5, 10, 15, 20, Inf)
  bmin_exp <- 1
  lambda_exp <- 0.1
  
  prden_exp <- getPDF(boundaries_exp, "expn", bmin_exp, lambda_exp)
  validate_pdf_output(prden_exp, boundaries_exp, "Exponential")
  
  # Test formula specifica per exponential
  # P(bi <= X < bi+1 | X >= bmin) = [exp(-λ*bi) - exp(-λ*bi+1)] / exp(-λ*bmin)
  l_exp <- boundaries_exp[-length(boundaries_exp)]
  u_exp <- boundaries_exp[-1]
  u_exp[is.infinite(u_exp)] <- 100  # Approssimazione per test
  
  theoretical_exp <- (exp(-lambda_exp * l_exp) - exp(-lambda_exp * u_exp)) / exp(-lambda_exp * bmin_exp)
  theoretical_exp <- theoretical_exp / sum(theoretical_exp)
  
  rel_errors_exp <- abs(prden_exp - theoretical_exp) / theoretical_exp
  if (any(rel_errors_exp > 0.01, na.rm = TRUE)) {
    warning("Exponential: Large relative errors in formula verification")
  }
  
  # ===================================================================
  # TEST 3: LOG-NORMAL
  # ===================================================================
  cat("Test 3: Log-normal Distribution\n")
  
  boundaries_ln <- c(0.1, 1, 5, 10, 50, 100)
  bmin_ln <- 0.1
  mu_ln <- 1.0
  sigma_ln <- 0.5
  
  prden_ln <- getPDF(boundaries_ln, "lgnorm", bmin_ln, mu_ln, sigma_ln)
  validate_pdf_output(prden_ln, boundaries_ln, "Log-normal")
  
  # Controllo: sigma > 0 
  if (sigma_ln <= 0) {
    stop("Log-normal test setup error: sigma must be > 0")
  }
  
  # ===================================================================
  # TEST 4: STRETCHED EXPONENTIAL (WEIBULL)
  # ===================================================================
  cat("Test 4: Stretched Exponential Distribution\n")
  
  boundaries_se <- c(1, 2, 5, 10, 20, 50)
  bmin_se <- 1
  lambda_se <- 0.1
  beta_se <- 0.8
  
  prden_se <- getPDF(boundaries_se, "stexp", bmin_se, lambda_se, beta_se)
  validate_pdf_output(prden_se, boundaries_se, "Stretched-exponential")
  
  # Controllo: parametri positivi
  if (lambda_se <= 0 || beta_se <= 0) {
    stop("Stretched exponential test setup error: lambda > 0 and beta > 0 required")
  }
  
  # ===================================================================
  # TEST 5: POWER-LAW CON CUTOFF ESPONENZIALE
  # ===================================================================
  cat("Test 5: Power-law with Exponential Cutoff\n")
  
  boundaries_plc <- c(1, 5, 10, 20, 50, 100)
  bmin_plc <- 1
  alpha_plc <- 2.5
  lambda_plc <- 0.05
  
  prden_plc <- getPDF(boundaries_plc, "plcut", bmin_plc, alpha_plc, lambda_plc)
  validate_pdf_output(prden_plc, boundaries_plc, "Power-law+cutoff")
  
  # Controllo: parametri validi
  if (alpha_plc <= 1 || lambda_plc <= 0) {
    stop("Power-law+cutoff test setup error: alpha > 1 and lambda > 0 required")
  }
  
  # ===================================================================
  # TEST 6: EDGE CASES E ROBUSTNESS
  # ===================================================================
  cat("Test 6: Edge Cases and Robustness\n")
  
  # Test con pochi bin
  boundaries_small <- c(1, 10, 100) # due bin
  prden_small <- getPDF(boundaries_small, "pl", 1, 2.0) # alpha = 2
  validate_pdf_output(prden_small, boundaries_small, "Small-boundaries")
  
  # Test con bin molto ampi
  boundaries_wide <- c(1, 1000, 1000000)
  prden_wide <- getPDF(boundaries_wide, "expn", 1, 0.001) # λ = 0.001
  validate_pdf_output(prden_wide, boundaries_wide, "Wide-boundaries")
  
  # Test gestione errori 
  tryCatch({
    # Alpha <= 1 per power-law
    invalid_result <- getPDF(c(1, 10, 100), "pl", 1, 0.5)
    if (!is.null(invalid_result)) {
      warning("getPDF should handle alpha <= 1 case")
    }
  }, error = function(e) {
    cat("Correctly caught error for invalid alpha\n")
  })
  
  tryCatch({
    # Sigma <= 0 per log-normal
    invalid_result <- getPDF(c(1, 10, 100), "lgnorm", 1, 0, -0.5)
    if (!is.null(invalid_result)) {
      warning("getPDF should handle sigma <= 0 case")
    }
  }, error = function(e) {
    cat("Correctly caught error for invalid sigma\n")
  })
  
  # ===================================================================
  # TEST 7: CONSISTENCY ACROSS DIFFERENT BMIN
  # ===================================================================
  cat("Test 7: Consistency with Different bmin Values\n")
  
  boundaries_cons <- c(1, 5, 10, 20, 50, 100)
  
  # Test con bmin = 1 (usa tutti i boundaries)
  prden_bmin1 <- getPDF(boundaries_cons, "pl", 1, 2.0)
  validate_pdf_output(prden_bmin1, boundaries_cons, "Consistency-bmin1")
  
  # Test con bmin = 5 (usa solo i boundaries >= 5)
  bmin_5 <- 5
  boundaries_for_bmin5 <- boundaries_cons[boundaries_cons >= bmin_5] # Filtra i boundaries
  prden_bmin5 <- getPDF(boundaries_for_bmin5, "pl", bmin_5, 2.0)
  
  validate_pdf_output(prden_bmin5, boundaries_for_bmin5, "Consistency-bmin5") # Usa i boundaries filtrati per la validazione
  
  return(TRUE)
}

# ===================================================================
# FUNZIONE DI TEST PARAMETRICA
# ===================================================================

test_getPDF_parametric <- function() {
  cat("=== PARAMETRIC TEST getPDF ===\n")
  
  # Test cases as list
  test_cases <- list(
    list(name = "PowerLaw_Basic", boundaries = c(1, 10, 100, 1000), 
         type = "pl", bmin = 1, params = list(2.0)),
    
    list(name = "PowerLaw_Steep", boundaries = c(1, 2, 4, 8, 16), 
         type = "pl", bmin = 1, params = list(3.5)),
    
    list(name = "Exponential_Fast", boundaries = c(0, 5, 10, 15, 20), 
         type = "expn", bmin = 0, params = list(0.5)),
    
    list(name = "Exponential_Slow", boundaries = c(1, 10, 20, 30, 40), 
         type = "expn", bmin = 1, params = list(0.05)),
    
    list(name = "LogNormal_Narrow", boundaries = c(0.1, 1, 5, 10, 20), 
         type = "lgnorm", bmin = 0.1, params = list(1.0, 0.3)),
    
    list(name = "LogNormal_Wide", boundaries = c(0.1, 1, 10, 100, 1000), 
         type = "lgnorm", bmin = 0.1, params = list(2.0, 1.5)),
    
    list(name = "StretchedExp", boundaries = c(1, 3, 6, 10, 15, 25), 
         type = "stexp", bmin = 1, params = list(0.1, 0.7)),
    
    list(name = "PowerLawCutoff", boundaries = c(1, 5, 10, 25, 50, 100), 
         type = "plcut", bmin = 1, params = list(2.2, 0.02))
  )
  
  # Run all test cases
  for (test_case in test_cases) {
    cat(sprintf("Testing: %s\n", test_case$name))
    
    # Call getPDF with appropriate parameters
    if (length(test_case$params) == 1) {
      result <- getPDF(test_case$boundaries, test_case$type, test_case$bmin, test_case$params[[1]])
    } else if (length(test_case$params) == 2) {
      result <- getPDF(test_case$boundaries, test_case$type, test_case$bmin, 
                       test_case$params[[1]], test_case$params[[2]])
    }
    
    # Validate
    if (is.null(result)) {
      stop(sprintf("Test %s failed: getPDF returned NULL", test_case$name))
    }
    
    expected_length <- length(test_case$boundaries) - 1
    if (length(result) != expected_length) {
      stop(sprintf("Test %s failed: wrong length", test_case$name))
    }
    
    if (any(result < 0)) {
      stop(sprintf("Test %s failed: negative probabilities", test_case$name))
    }
    
    sum_result <- sum(result)
    if (abs(sum_result - 1.0) > 1e-8) {
      stop(sprintf("Test %s failed: poor normalization (sum=%.10f)", test_case$name, sum_result))
    }
    
    cat(sprintf("%s: OK (sum=%.8f)\n", test_case$name, sum_result))
  }
  
  cat("All parametric tests passed\n")
  return(TRUE)
}

# ===================================================================
# ESECUZIONE DEI TEST
# ===================================================================

if (interactive()) {
  cat("Running getPDF comprehensive tests...\n\n")
  
  # Test principale
  test_getPDF_complete()
  
  # stampa ogni argomento in sequenza
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  
  # Test parametrico
  test_getPDF_parametric()
  
  cat("\nAll getPDF tests completed successfully!\n")
}