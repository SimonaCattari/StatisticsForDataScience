# ===============================================================================
# FRAMEWORK VIRKAR & CLAUSET (2014) - POWER LAWS IN BINNED DATA
# ===============================================================================

library(httr)
library(readr) # lettura file di testo

cat("FRAMEWORK VIRKAR & CLAUSET (2014)\n")
cat("===================================\n")
cat("Power-law analysis su dati binned\n")
cat("Dataset: Forest Fires\n\n")

# ===============================================================================
# 1. CARICAMENTO FUNZIONI
# ===============================================================================

cat("Caricamento funzioni del framework...\n")

base_path <- "Funzioni" 

# Lista funzioni necessarie
required_functions <- c(
  "bplfit.R",      # Fit power law
  "getPDF.R",      # Calcolo PDF
  "bplpva.R",      # Test plausibilità
  "blrtest.R",     # Likelihood ratio test
  "bexpnfit.R",    # Fit exponential
  "bstexpfit.R",   # Fit stretched exponential
  "blgnormfit.R",  # Fit log-normal
  "bplcutfit.R"    # Fit power law + cutoff
)

# Carica funzioni
functions_loaded <- 0
for (func in required_functions) {
  func_path <- file.path(base_path, func)
  if (file.exists(func_path)) {
    tryCatch({
      source(func_path)
      cat(sprintf("Caricate%s\n", func))
      functions_loaded <- functions_loaded + 1
    }, error = function(e) {
      cat(sprintf("Errors: %s\n", func, e$message))
    })
  } else {
    cat(sprintf("%s: file non trovato\n", func))
  }
}

if (functions_loaded < 6) {
  stop("Errore: troppe poche funzioni caricate. Verifica il path.")
}

cat(sprintf("%d/%d funzioni caricate con successo\n\n", functions_loaded, length(required_functions)))

# ===============================================================================
# 2. FUNZIONI HELPER
# ===============================================================================

# Funzione per fit power law
fit_power_law <- function(h, boundaries) {
  tryCatch({
    bplfit(h, boundaries)
  }, error = function(e) {
    cat(sprintf("Errore in bplfit: %s\n", e$message))
    list(alpha = NA, bmin = NA, logLik = NA, D = NA)
  })
}

# Funzione per fit alternative
fit_alternatives <- function(h, boundaries, bmin_pl, alpha_pl) {
  
  cat("Fitting alternative (con range generici)...\n")
  
  results <- list()
  
  # Prepara dati per fit (solo sulla coda)
  # Trova l'indice del bmin all'interno delle boundaries
  bmin_idx <- max(which(boundaries <= bmin_pl))
  if (is.infinite(bmin_idx) || bmin_idx < 1 || bmin_idx >= length(boundaries)) {
    warning("bmin_pl non valido. Usando tutti i dati per fit alternativi.")
    h_tail <- h
    boundaries_tail <- boundaries
  } else {
    h_tail <- h[bmin_idx:length(h)]
    boundaries_tail <- boundaries[bmin_idx:length(boundaries)]
  }
  
  # 1. EXPONENTIAL
  cat("Exponential...\n")
  results$exp <- tryCatch({
    # NOTA: Range di lambda molto ampio per adattarsi a scale diverse
    lambda_range <- c(10^seq(-10, 0, length.out=100))
    result <- bexpnfit(h_tail, boundaries_tail, lambda_range = lambda_range, bmin = bmin_pl)
    if (!is.na(result$lambda) && result$lambda > 0 && is.finite(result$loglik)) {
      cat(sprintf("λ = %.8f, logLik = %.2f\n", result$lambda, result$loglik))
      result
    } else { list(lambda = NA, loglik = NA) }
  }, error = function(e) {
    cat(sprintf("Errore in Exponential fit: %s\n", e$message))
    list(lambda = NA, loglik = NA)
  })
  
  # 2. LOG-NORMAL
  cat("Log-normal...\n")
  results$lognorm <- tryCatch({
    mu_range <- seq(0, 15, by = 0.2)
    sigma_range <- seq(0.1, 5.0, by = 0.1)
    result <- blgnormfit(h_tail, boundaries_tail, murng = mu_range, sigrng = sigma_range, bmin = bmin_pl, fine = TRUE)
    if (!is.na(result$mu) && !is.na(result$sigma) && result$sigma > 0 && is.finite(result$logLik)) {
      cat(sprintf("μ = %.3f, σ = %.3f, logLik = %.2f\n", result$mu, result$sigma, result$logLik))
      result
    } else { list(mu = NA, sigma = NA, logLik = NA, bmin = bmin_pl) }
  }, error = function(e) {
    cat(sprintf("Errore in Log-normal fit: %s\n", e$message))
    list(mu = NA, sigma = NA, logLik = NA, bmin = bmin_pl)
  })
  
  # 3. POWER LAW + CUTOFF
  cat("Power law + cutoff...\n")
  results$plcut <- tryCatch({
    alpha_range <- seq(1.2, 4.0, by = 0.1)
    lambda_range <- c(10^seq(-12, -1, length.out = 30))
    result <- bplcutfit(h_tail, boundaries_tail, range_alpha = alpha_range, range_lambda = lambda_range, bmin = bmin_pl)
    if (!is.na(result$alpha) && !is.na(result$lambda) && result$alpha > 1 && result$lambda > 0 && is.finite(result$loglik)) {
      cat(sprintf("α = %.3f, λ = %.10f, logLik = %.2f\n", result$alpha, result$lambda, result$loglik))
      result
    } else { list(alpha = NA, lambda = NA, loglik = NA) }
  }, error = function(e) {
    cat(sprintf("Errore in PL+cutoff fit: %s\n", e$message))
    list(alpha = NA, lambda = NA, loglik = NA)
  })
  
  # 4. STRETCHED EXPONENTIAL
  cat("Stretched exponential...\n")
  results$strexp <- tryCatch({
    lambda_range <- c(10^seq(-15, 0, length.out = 40))
    beta_range <- seq(0.1, 2.5, by = 0.1)
    result <- bstexpfit(h_tail, boundaries_tail, range = list(lambda = lambda_range, beta = beta_range), bmin = bmin_pl)
    if (!is.na(result$lambda) && !is.na(result$beta) && result$lambda > 0 && result$beta > 0 && is.finite(result$loglik)) {
      cat(sprintf("λ = %.12f, β = %.3f, logLik = %.2f\n", result$lambda, result$beta, result$loglik))
      result
    } else { list(lambda = NA, beta = NA, loglik = NA) }
  }, error = function(e) {
    cat(sprintf("Errore in Stretched Exponential fit: %s\n", e$message))
    list(lambda = NA, beta = NA, loglik = NA)
  })
  
  return(results)
}

# Stampa i risultati LRT
print_lr_results <- function(lr_results) {
  cat("\nRISULTATI LIKELIHOOD RATIO TESTS:\n")
  
  # vs Log-normal
  if (!is.null(lr_results$lognorm) && !is.na(lr_results$lognorm$normR)) {
    cat(sprintf("vs Log-normal: R = %.3f, p = %.3f", lr_results$lognorm$normR, lr_results$lognorm$p))
    if (!is.na(lr_results$lognorm$p) && lr_results$lognorm$p < 0.1) {
      if (lr_results$lognorm$normR > 0) cat(" -> PL favorita") else cat(" -> Log-normal favorita")
    } else cat(" -> Indeciso")
    cat("\n")
  } else { cat("vs Log-normal: FALLITO\n") }
  
  # vs Exponential
  if (!is.null(lr_results$exp) && !is.na(lr_results$exp$normR)) {
    cat(sprintf("vs Exponential: R = %.3f, p = %.3f", lr_results$exp$normR, lr_results$exp$p))
    if (!is.na(lr_results$exp$p) && lr_results$exp$p < 0.1) {
      if (lr_results$exp$normR > 0) cat(" -> PL favorita") else cat(" -> Exponential favorita")
    } else cat(" -> Indeciso")
    cat("\n")
  } else { cat("vs Exponential: FALLITO\n") }
  
  # vs Stretched exponential
  if (!is.null(lr_results$strexp) && !is.na(lr_results$strexp$normR)) {
    cat(sprintf("vs Stretched exp: R = %.3f, p = %.3f", lr_results$strexp$normR, lr_results$strexp$p))
    if (!is.na(lr_results$strexp$p) && lr_results$strexp$p < 0.1) {
      if (lr_results$strexp$normR > 0) cat(" -> PL favorita") else cat(" -> Stretched exp favorita")
    } else cat(" -> Indeciso")
    cat("\n")
  } else { cat("vs Stretched exp: FALLITO\n") }
  
  # vs Power law + cutoff
  if (!is.null(lr_results$plcut) && !is.na(lr_results$plcut$normR)) {
    cat(sprintf("vs PL+cutoff: R = %.3f, p = %.3f", lr_results$plcut$normR, lr_results$plcut$p))
    if (!is.na(lr_results$plcut$p) && lr_results$plcut$p < 0.1) {
      if (lr_results$plcut$normR > 0) cat(" -> PL favorita") else cat(" -> PL+cutoff favorita")
    } else cat(" -> Indeciso")
    cat("\n")
  } else { cat("vs PL+cutoff: DFALLITO\n") }
}

# Funzione il livello di supporto
determine_support_level <- function(fit_success, plausible, lr_results) {
  if (!fit_success) return("none")
  if (!plausible) return("none") 
  
  strong_alternatives <- 0
  
  # Controllo  se le alternative sono significativamente migliori (p < 0.1 e R < 0)
  if (!is.null(lr_results$lognorm) && !is.na(lr_results$lognorm$p) && lr_results$lognorm$p < 0.1 && lr_results$lognorm$normR < 0) {
    strong_alternatives <- strong_alternatives + 1
  }
  if (!is.null(lr_results$strexp) && !is.na(lr_results$strexp$p) && lr_results$strexp$p < 0.1 && lr_results$strexp$normR < 0) {
    strong_alternatives <- strong_alternatives + 1
  }
  
  # PL+Cutoff è un caso speciale (modello nested)
  is_plcut_better <- !is.null(lr_results$plcut) && !is.na(lr_results$plcut$p) && lr_results$plcut$p < 0.1 && lr_results$plcut$normR < 0
  
  if (strong_alternatives > 0) {
    return("weak")
  }
  if (is_plcut_better) {
    return("moderate") # Se PL+Cutoff è l'unica alternativa migliore, il supporto è moderato
  }
  
  return("good") # Se la PL è plausibile e nessuna alternativa è significativamente migliore
}

# Funzione per eseguire i test di likelihood
perform_lr_tests <- function(h, boundaries, bmin, alpha_pl, alternative_fits) {
  lr_results <- list()
  
  pdf_pl <- tryCatch(getPDF(boundaries, "pl", bmin, alpha_pl), error = function(e) NULL)
  if (is.null(pdf_pl)) return(lr_results)
  
  # vs Log-normal
  if (!is.null(alternative_fits$lognorm) && !is.na(alternative_fits$lognorm$mu)) {
    pdf_ln <- tryCatch(getPDF(boundaries, "lgnorm", bmin, alternative_fits$lognorm$mu, alternative_fits$lognorm$sigma), error = function(e) NULL)
    if (!is.null(pdf_ln)) lr_results$lognorm <- tryCatch(blrtest(pdf_pl, pdf_ln, h, boundaries, bmin), error = function(e) list(normR = NA, p = NA))
  }
  
  # vs Exponential
  if (!is.null(alternative_fits$exp) && !is.na(alternative_fits$exp$lambda)) {
    pdf_exp <- tryCatch(getPDF(boundaries, "expn", bmin, alternative_fits$exp$lambda), error = function(e) NULL)
    if (!is.null(pdf_exp)) lr_results$exp <- tryCatch(blrtest(pdf_pl, pdf_exp, h, boundaries, bmin), error = function(e) list(normR = NA, p = NA))
  }
  
  # vs Stretched exponential
  if (!is.null(alternative_fits$strexp) && !is.na(alternative_fits$strexp$lambda)) {
    pdf_strexp <- tryCatch(getPDF(boundaries, "stexp", bmin, alternative_fits$strexp$lambda, alternative_fits$strexp$beta), error = function(e) NULL)
    if (!is.null(pdf_strexp)) lr_results$strexp <- tryCatch(blrtest(pdf_pl, pdf_strexp, h, boundaries, bmin), error = function(e) list(normR = NA, p = NA))
  }
  
  # vs Power law + cutoff
  if (!is.null(alternative_fits$plcut) && !is.na(alternative_fits$plcut$alpha)) {
    pdf_plcut <- tryCatch(getPDF(boundaries, "plcut", bmin, alternative_fits$plcut$alpha, alternative_fits$plcut$lambda), error = function(e) NULL)
    if (!is.null(pdf_plcut)) lr_results$plcut <- tryCatch(blrtest(pdf_pl, pdf_plcut, h, boundaries, bmin, isNested = TRUE), error = function(e) list(normR = NA, p = NA))
  }
  
  return(lr_results)
}

# ===============================================================================
# 3. DOWNLOAD E PREPARAZIONE DATASET
# ===============================================================================

cat("\nDownload dataset 'Forest Fires'...\n")

dataset_url <- "https://sites.santafe.edu/~aaronc/powerlaws/data/fires.txt"

fire_data <- tryCatch({
  temp_file <- tempfile(fileext = ".txt")
  download.file(dataset_url, temp_file, mode = "wb", quiet = TRUE)
  raw_data <- readLines(temp_file, warn = FALSE)
  unlink(temp_file)
  
  # Converti a numerico e pulisci
  values <- as.numeric(raw_data)
  values <- values[!is.na(values) & values > 0]
  
  cat(sprintf("%d valori scaricati dal dataset\n", length(values)))
  sort(values, decreasing = TRUE)
  
}, error = function(e) {
  stop(sprintf("Download fallito: %s\n", e$message))
})

cat(sprintf("Dataset processato:\n"))
cat(sprintf(" -> N valori: %s\n", format(length(fire_data), big.mark = ",")))
cat(sprintf(" -> Range: %s - %s\n",
            format(min(fire_data), big.mark = ","),
            format(max(fire_data), big.mark = ",")))

# ===============================================================================
# 4. BINNING LOGARITMICO
# ===============================================================================

cat("\nBinning logaritmico (c = 2)...\n")

min_val <- min(fire_data)
max_val <- max(fire_data)

start_power <- floor(log2(min_val))
end_power <- ceiling(log2(max_val)) + 1
boundaries <- 2^(start_power:end_power)

h <- hist(fire_data, breaks = boundaries, plot = FALSE, right = FALSE)$counts

cat(sprintf("Binning completato:\n"))
cat(sprintf(" -> Boundaries: %d (da 2^%d a 2^%d)\n", length(boundaries), start_power, end_power-1))
cat(sprintf(" -> Bins non vuoti: %d\n", sum(h > 0)))
cat(sprintf(" -> Verifica: %d osservazioni totali\n", sum(h)))

# ===============================================================================
# 5. STEP 1: FIT POWER LAW
# ===============================================================================

cat("\nSTEP 1: FIT POWER LAW\n")
cat("========================\n")

pl_fit <- fit_power_law(h, boundaries)

if (!is.na(pl_fit$alpha) && pl_fit$alpha > 1 && !is.na(pl_fit$bmin) && pl_fit$bmin > 0) {
  bmin_idx_for_tail <- max(which(boundaries <= pl_fit$bmin))
  h_tail <- if (bmin_idx_for_tail <= length(h)) h[bmin_idx_for_tail:length(h)] else 0
  n_tail <- sum(h_tail)
  
  cat("RISULTATI STEP 1:\n")
  cat(sprintf(" -> α = %.3f\n", pl_fit$alpha))
  cat(sprintf(" -> bmin = %s\n", format(pl_fit$bmin, big.mark = ",")))
  cat(sprintf(" -> Log-likelihood = %.2f\n", pl_fit$logLik))
  cat(sprintf(" -> KS statistic = %.6f\n", pl_fit$D))
  cat(sprintf(" -> Osservazioni coda (n_tail) = %d (%.1f%%)\n", n_tail, 100 * n_tail / sum(h)))
  fit_success <- TRUE
} else {
  cat("STEP 1 FALLITO: parametri non validi\n")
  fit_success <- FALSE
  n_tail <- NA
}

# ===============================================================================
# 6. STEP 2: TEST PLAUSIBILITÀ
# ===============================================================================

cat("\n STEP 2: TEST PLAUSIBILITÀ\n")
cat("============================\n")

if (fit_success) {
  cat("Eecuzione bootstrap test (1000 repliche)..\n")
  pva_result <- tryCatch({
    bplpva(h, boundaries, pl_fit$bmin, pl_fit$alpha, reps = 1000, silent = TRUE, seed = 123)
  }, error = function(e) {
    cat(sprintf("plpva fallito: %s\n", e$message))
    list(p = NA, Dstar = NA, n_valid_bootstrap = 0, d_vals = NA)
  })
  
  if (!is.na(pva_result$p)) {
    cat("RISULTATI STEP 2:\n")
    cat(sprintf(" -> P-value = %.3f\n", pva_result$p))
    if (pva_result$p >= 0.1) {
      cat("POWER LAW PLAUSIBILE (p ≥ 0.1)\n")
      plausible <- TRUE
    } else {
      cat("POWER LAW RIGETTATA (p < 0.1)\n")
      plausible <- FALSE
    }
  } else {
    cat("STEP 2 FALLITO\n")
    plausible <- FALSE
  }
} else {
  cat("STEP 2 SALTATO\n")
  plausible <- FALSE
  pva_result <- list(p = NA) 
}

# ===============================================================================
# 7. STEP 3: CONFRONTO ALTERNATIVE
# ===============================================================================

cat("\nSTEP 3: CONFRONTO ALTERNATIVE\n")
cat("=================================\n")

alternative_fits <- list()
lr_results <- list()

if (fit_success) {
  alternative_fits <- fit_alternatives(h, boundaries, pl_fit$bmin, pl_fit$alpha)
  lr_results <- perform_lr_tests(h, boundaries, pl_fit$bmin, pl_fit$alpha, alternative_fits)
  print_lr_results(lr_results)
} else {
  cat("STEP 3 SALTATO\n")
}

# ===============================================================================
# 8. VALUTAZIONE FINALE
# ===============================================================================

cat("\nVALUTAZIONE FINALE\n")
cat("=====================\n")

support_level <- determine_support_level(fit_success, plausible, lr_results)
cat(sprintf("SUPPORTO PER POWER LAW: %s\n", toupper(support_level)))

# ===============================================================================
# 9. VISUALIZZAZIONI
# ===============================================================================

cat("\nCreazione visualizzazioni...\n")

tryCatch({
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))
  
  bin_centers <- sqrt(boundaries[-length(boundaries)] * boundaries[-1])
  empirical_ccdf <- rev(cumsum(rev(h))) / sum(h)
  valid_idx <- empirical_ccdf > 0 & bin_centers > 0
  
  # 1. CCDF plot
  plot(log10(bin_centers[valid_idx]), log10(empirical_ccdf[valid_idx]),
       pch = 16, col = "darkblue", cex = 0.8,
       xlab = "log10(Area Incendio)", ylab = "log10(P(X ≥ x))",
       main = "CCDF Empirica vs Power Law")
  
  if (fit_success) {
    x_fit <- seq(pl_fit$bmin, max(boundaries), length.out = 100)
    bmin_ccdf_val <- empirical_ccdf[which.min(abs(bin_centers - pl_fit$bmin))]
    y_fit_ccdf <- (x_fit / pl_fit$bmin)^(-(pl_fit$alpha - 1)) * bmin_ccdf_val
    lines(log10(x_fit), log10(y_fit_ccdf), col = "red", lwd = 2)
    abline(v = log10(pl_fit$bmin), col = "red", lty = 2)
    legend("bottomleft", c("Dati", "Power law", "bmin"), col = c("darkblue", "red", "red"), pch = c(16, NA, NA), lty = c(NA, 1, 2), cex = 0.7)
  }
  
  # 2. Istogramma binned
  plot(log10(bin_centers), h, type = "h", lwd = 2, col = "steelblue",
       xlab = "log10(Area Incendio)", ylab = "Frequenza",
       main = "Istogramma Binned e Fit")
  if (fit_success) {
    abline(v = log10(pl_fit$bmin), col = "red", lty = 2, lwd = 2)
    legend("topright", c("Dati binned", "bmin"), col = c("steelblue", "red"), lty = c(1, 2), pch = c(15, NA), cex = 0.7)
  }
  
  # 3. Distribuzione bootstrap
  if (exists("pva_result") && !is.null(pva_result$d_vals)) {
    hist(pva_result$d_vals, breaks = 30, col = "lightsteelblue", main = "Bootstrap Distribution", xlab = "D (KS statistic)")
    abline(v = pva_result$Dstar, col = "red", lwd = 3)
    title(sub = paste("p-value =", round(pva_result$p, 3)))
  } else {
    plot(1, 1, type = "n", main = "Bootstrap Non Disponibile", xlab = "", ylab = ""); text(1, 1, "Test non eseguito", col = "red")
  }
  
  # 4. Summary risultati
  plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  summary_text <- c(
    paste("Dataset:", length(fire_data), "incendi"), "", "RISULTATI:",
    if(fit_success) paste("α =", round(pl_fit$alpha, 3)) else "α = FALLITO",
    if(fit_success) paste("bmin =", format(round(pl_fit$bmin, 1), big.mark = ",")) else "bmin = FALLITO",
    if(!is.na(pva_result$p)) paste("p-value =", round(pva_result$p, 3)) else "p-value = FALLITO", "",
    paste("SUPPORTO:", toupper(support_level))
  )
  text(0.05, 0.9, paste(summary_text, collapse = "\n"), adj = c(0, 1), cex = 0.9, family = "mono")
  title("Riepilogo Analisi: Forest Fires", outer = TRUE)
  par(mfrow = c(1, 1))
  cat("Visualizzazioni create\n")
  
}, error = function(e) {
  cat(sprintf("Errore visualizzazioni: %s\n", e$message))
})

# ===============================================================================
# 10. CONCLUSIONE
# ===============================================================================

cat("\n FRAMEWORK COMPLETATO\n")
cat("=======================\n")
if(fit_success) {
  cat(sprintf(" -> Power Law: α=%.3f, bmin=%s\n", pl_fit$alpha, format(round(pl_fit$bmin,1), big.mark = ",")))
  cat(sprintf(" -> Plausibilità: %s (p=%.3f)\n", ifelse(plausible, "SÌ", "NO"), ifelse(is.na(pva_result$p), 0, pva_result$p)))
  cat(sprintf(" -> Supporto Statistico Finale: %s\n", toupper(support_level)))
} else {
  cat(" -> Analisi FALLITA per problemi nel fit iniziale.\n")
}
cat("\nAnalisi del dataset 'fires.txt' completata.\n")