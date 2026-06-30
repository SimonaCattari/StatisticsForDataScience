# =============================================================================
# - Stima il parametro di scala alpha tramite massimo di verosimiglianza (MLE)
# - Identifica automaticamente il valore ottimale di bmin, cioè il punto a partire
#   dal quale si assume che la distribuzione segua una power-law
# - Utilizza la statistica di Kolmogorov-Smirnov (KS) per selezionare il miglior bmin
#
# La funzione può essere anche vincolata a un range di valori di alpha (`range`),
# a un limite massimo per bmin (`limit`), oppure a un valore fissato di bmin (`bmin_fixed`).
#
# ARGOMENTI:
#   - h:              Vettore di conteggi per ciascun bin (numeri interi non negativi)
#   - boundaries:     Vettore degli estremi dei bin (lunghezza = length(h) + 1)
#   - range:          Intervallo di ricerca per il parametro alpha (es. c(1.5, 3)) (opzionale)
#   - limit:          Valore massimo consentito per bmin (opzionale)
#   - bmin_fixed:     Valore prefissato di bmin da utilizzare (opzionale)
#
# OUTPUT:
#   - alpha:    stima MLE del parametro di scala
#   - bmin:     valore di soglia per la coda (bin minimo)
#   - logLik:   log-verosimiglianza del modello power-law sopra bmin
#   - D:        valore della statistica KS tra dati osservati e modello teorico
#
# =============================================================================

bplfit <- function(h, boundaries, range = NULL, limit = NULL, bmin_fixed = NULL) {
  if (!all(h == floor(h)) || any(h < 0)) stop("'h' must be a non negtive integer.")
  if (length(boundaries) != length(h) + 1) stop("'boundaries' must have length(h) + 1")
  if (length(h) < 2) stop("At least 2 bins required")
  
  # Si mantengono almeno 2 bin significativi e si rimuovono 
  # i bin vuoti per migliorare la stabilità
  while (length(h) > 2 && h[length(h)] == 0) {
    h <- h[-length(h)]
    boundaries <- boundaries[-(length(boundaries))]
  }
  
  # Determinare i candidati per bmin. Servono almeno 2 bin dopo b_min, 
  # quindi escludere gli ultimi due estremi.
  bmin_candidates <- boundaries[1:(length(boundaries)-2)] 
  
  # Si restringono i candidati per bmin in base ai limiti impostati (se esistono)
  if (!is.null(limit)) {
    bmin_candidates <- bmin_candidates[bmin_candidates <= limit]
  }
  if (!is.null(bmin_fixed)) {
    # Se viene utilizzato bmin_fixed, testiamo solo i candidati uguali o vicini a bmin_fixed.
    # Per semplicità, qui consideriamo solo bmin_fixed se coincide con uno degli estremi (boundaries).
    if (!(bmin_fixed %in% boundaries)) {
      warning("'bmin_fixed' is not one of the bin boundaries. Using closest boundary.")
      closest_bmin_candidate <- boundaries[which.min(abs(boundaries - bmin_fixed))]
      bmin_candidates <- bmin_candidates[bmin_candidates == closest_bmin_candidate]
    } else {
      bmin_candidates <- bmin_candidates[bmin_candidates == bmin_fixed]
    }
  }
  
  if (length(bmin_candidates) == 0) {
    return(list(alpha = NA, bmin = NA, logLik = -Inf, D = NA))
  }
  
  alphas <- rep(NA, length(bmin_candidates)) # stima di α per ogni candidato bmin
  logLiks <- rep(NA, length(bmin_candidates)) # valori della log-verosimiglianza calcolata via MLE
  ks_stats <- rep(NA, length(bmin_candidates)) # valori della statistica di K-S
  
  # MLE per stimare α
  compute_alpha_mle <- function(h_tail, l_bounds, u_bounds, bmin_current) {
    n_tail <- sum(h_tail)
    if (n_tail == 0) return(list(alpha = NA, logLik = -Inf))
    
    neg_loglik <- function(alpha) { # definita solo per α > 1
      if (alpha <= 1) return(Inf)
      
      # Membri equazione 3.1
      prob_unnorm_num <- l_bounds^(1-alpha) - u_bounds^(1-alpha) # Numeratore
      norm_denom <- bmin_current^(1-alpha) # Denominatore

      log_bin_probs <- log(pmax(prob_unnorm_num / norm_denom, .Machine$double.xmin)) # Usa un piccolo epsilon per evitare log(0)
      
      # Se qualche bin con h > 0 ha log_prob non finito (NaN o -Inf), restituisce Inf: α non valido
      if (any(h_tail[h_tail > 0] > 0 & (!is.finite(log_bin_probs[h_tail > 0])))) return(Inf)
      
      # Formula completa dall'Eq. 3.1 
      loglik <- sum(h_tail * log_bin_probs)
      
      return(-loglik)
    }
    
    # Ricerca del miglior paramentro α:
    # 1) Grid Search
    if (!is.null(range)) {
      alpha_grid <- seq(min(range), max(range), length.out = 100) # Si provano 100 valori equidistanti tra i limiti
      nll_values <- sapply(alpha_grid, neg_loglik)
      
      if (all(is.infinite(nll_values))) return(list(alpha = NA, logLik = -Inf))
      
      best_alpha_idx <- which.min(nll_values)
      best_alpha <- alpha_grid[best_alpha_idx]
      best_loglik <- -nll_values[best_alpha_idx]
      # Prende il valore di α che minimizza la neg_loglik
      
      return(list(alpha = best_alpha, logLik = best_loglik))
    } else {
      # 2) Ricerca più precisa con 'optimize()'
      opt_result <- tryCatch({
        optimize(neg_loglik, interval = c(1.0001, 100), tol = 1e-10) # intervallo α nell’intervallo [1.0001, 100]
      }, error = function(e) {
        warning(paste("Ottimizzazione fallita:", e$message))
        list(minimum = NA, objective = Inf)
      })
      
      if (is.finite(opt_result$minimum) && opt_result$minimum > 1) {
        return(list(alpha = opt_result$minimum, logLik = -opt_result$objective))
      } else {
        return(list(alpha = NA, logLik = -Inf))
      }
    }
  }
  
  # Iterazioni attraverso ciascun candidato bmin per trovare il miglior fit
  for (i in seq_along(bmin_candidates)) {
    bmin_candidate <- bmin_candidates[i]
    
    # Trova l'indice del primo bin il cui limite inferiore è <= bmin_candidate,
    # considera anche il caso in cui bmin_candidate sia tra due limiti
    bin_idx_start <- which(boundaries <= bmin_candidate)
    if (length(bin_idx_start) == 0) next
    bin_idx_start <- max(bin_idx_start) # limite massimo <= bmin_candidate
    
    # bmin_candidate è oltre l'ultimo bin valido
    if (bin_idx_start >= (length(boundaries) - 1)) next 
    
    h_tail <- h[bin_idx_start:(length(h))] # seleziona i dati nella coda dove si ipotizza la pl
    l_bounds_tail <- boundaries[bin_idx_start:(length(boundaries)-1)]
    u_bounds_tail <- boundaries[(bin_idx_start+1):length(boundaries)]
    
    # se non ci sono dati nella coda
    if (sum(h_tail) == 0) next 
    
    # Stima α con MLE per il 'bmin_candidate'
    alpha_mle_result <- compute_alpha_mle(h_tail, l_bounds_tail, u_bounds_tail, bmin_candidate)
    alpha_est <- alpha_mle_result$alpha
    
    if (is.na(alpha_est) || alpha_est <= 1) next
    
    alphas[i] <- alpha_est
    logLiks[i] <- alpha_mle_result$logLik 
    
    # K-S per il 'bmin_candidate
    n_tail <- sum(h_tail)
    ccdf_empirical <- cumsum(h_tail) / n_tail 
    ccdf_theoretical <- 1 - (u_bounds_tail / bmin_candidate)^(1 - alpha_est) # CDF teorica P(X < x) = 1 - (x/bmin)^(1-alpha)
    ccdf_theoretical <- pmax(0, pmin(1, ccdf_theoretical)) # Limita valori per evitare problemi numerivi
    
    # KS: massima differenza tra CDF empirica e teorica
    if (length(ccdf_empirical) != length(ccdf_theoretical)) {
      warning(paste("Length mismatch for KS stat at bmin_candidate:", bmin_candidate))
      ks_stats[i] <- NA
    } else {
      ks_stats[i] <- max(abs(ccdf_empirical - ccdf_theoretical))
    }
  }
  
  # Miglior 'bmin' (KS minima) 
  valid_indices <- which(!is.na(ks_stats) & is.finite(ks_stats))
  if (length(valid_indices) == 0) {
    return(list(alpha = NA, bmin = NA, logLik = -Inf, D = NA))
  }
  
  best_idx <- valid_indices[which.min(ks_stats[valid_indices])]
  best_bmin <- bmin_candidates[best_idx]
  best_alpha <- alphas[best_idx]
  best_ks <- ks_stats[best_idx]
  best_logLik_final <- logLiks[best_idx] 
  
  # Parametri migliori
  return(list(
    alpha = best_alpha,
    bmin = best_bmin,
    logLik = best_logLik_final, 
    D = best_ks
  ))
}