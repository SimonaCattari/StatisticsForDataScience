# =============================================================================
# bplpva(): Calcola il p-value per il test di bontà di adattamento (GoF)
#           di una distribuzione power-law su dati binned, tramite bootstrap
#           semi-parametrico.
#
# ARGOMENTI:
#   - h:           vettore di conteggi binned
#   - boundaries:  estremi dei bin (length = length(h) + 1)
#   - bmin_orig:   valore di bmin stimato (limite inferiore del power-law)
#   - alpha_orig:  esponente della power-law stimato (alpha > 1)
#   - reps:        numero di repliche bootstrap (default = 1000)
#   - silent:      se TRUE, sopprime output a console
#   - seed:        seme per la riproducibilità
#
# OUTPUT:
#   - p:           p-value del test KS
#   - d_vals:      vettore dei valori KS bootstrap
#   - Dstar:       statistica KS sui dati originali
#   - n_valid_bootstrap: numero di repliche bootstrap valide
#   - n_tail:      numero di osservazioni nella coda (tail)
# =============================================================================


source("Funzioni/bplfit.R") 

bplpva <- function(h, boundaries, bmin_orig, alpha_orig, reps = 1000, silent = FALSE, seed = NULL) {
  
  if (!is.numeric(h) || any(h < 0) || any(h != floor(h))) {
    stop("Errore: 'h' deve contenere conteggi interi non negativi.")
  }
  if (!is.numeric(boundaries) || length(boundaries) != length(h) + 1) {
    stop("Errore: 'boundaries' deve avere length(h) + 1 elementi.")
  }
  if (length(h) < 2) stop("Errore: 'h' deve contenere almeno 2 bin.")
  if (!is.numeric(bmin_orig) || length(bmin_orig) != 1 || bmin_orig <= 0) {
    stop("Errore: 'bmin_orig' deve essere uno scalare numerico positivo.")
  }
  if (!is.numeric(alpha_orig) || length(alpha_orig) != 1 || alpha_orig <= 1) {
    stop("Errore: 'alpha_orig' deve essere uno scalare numerico maggiore di 1.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # ===================================================================
  # SETUP PRINCIPALE
  # ===================================================================
  N <- sum(h)
  
  # Trova l'indice del primo bin >= bmin_orig
  idx_bmin <- which(boundaries[-length(boundaries)] >= bmin_orig)[1]
  if (is.na(idx_bmin)) stop("bmin_orig non trovato nei confini dei bin.")
  
  # Dati nella regione power-law (da bmin in poi)
  h_tail <- h[idx_bmin:length(h)]
  boundaries_tail <- boundaries[idx_bmin:length(boundaries)]
  
  if (sum(h_tail) < 2 || length(h_tail) < 2) {
    warning("Dati insufficienti nella regione power-law.")
    return(list(p = NA_real_, d_vals = rep(NA_real_, reps), Dstar = NA_real_))
  }
  
  n_tail <- sum(h_tail)
  
  # ===================================================================
  # CALCOLA D* ORIGINALE (KS SUI DATI EMPIRICI)
  # ===================================================================
  
  # CCDF empirica (Complementary CDF) per i dati nella coda
  # CCDF(x) = P(X >= x) = frazione di osservazioni >= x
  bin_lowers_tail <- boundaries_tail[-length(boundaries_tail)]  # Lower bounds dei bin
  
  # CCDF empirica: per ogni bin, conta quanti bin hanno valore >= quel bin
  ccdf_empirical <- rev(cumsum(rev(h_tail))) / n_tail
  
  # CCDF teorica per power law: P(X >= bi | X >= bmin) = (bi/bmin)^(1-α)
  ccdf_theoretical <- (bin_lowers_tail / bmin_orig)^(1 - alpha_orig)
  
  # Assicura che CCDF sia nel range [0,1]
  ccdf_theoretical <- pmax(0, pmin(1, ccdf_theoretical))
  
  # Statistica KS
  Dstar <- max(abs(ccdf_empirical - ccdf_theoretical))
  
  if (!silent) {
    cat(sprintf("D* originale = %.6f\n", Dstar))
    cat(sprintf("Range boundaries: [%.2f, %.2f]\n", boundaries[1], boundaries[length(boundaries)]))
  }
  
  # ===================================================================
  # BOOTSTRAP SEMI-PARAMETRICO
  # ===================================================================
  
  # Setup per bootstrap: dati sotto bmin
  if (idx_bmin > 1) {
    h_below <- h[1:(idx_bmin - 1)]
    below_boundaries <- boundaries[1:idx_bmin]
    # Centri dei bin sotto bmin
    mids_below <- (below_boundaries[-length(below_boundaries)] + below_boundaries[-1]) / 2
    x_pool_below <- rep(mids_below, times = h_below)
    N_below <- sum(h_below)
  } else {
    x_pool_below <- numeric(0)
    N_below <- 0
  }
  
  # Generazione power law continua
  rpl_continuous <- function(n, bmin, alpha) {
    if (n == 0) return(numeric(0))
    u <- runif(n)
    bmin * (1 - u)^(-1 / (alpha - 1))
  }
  
  # Binning 
  safe_binning <- function(data, boundaries) {
    if (length(data) == 0) {
      return(rep(0, length(boundaries) - 1))
    }
    
    # Estendi boundaries se necessario
    b_work <- boundaries
    if (min(data) < b_work[1]) {
      b_work[1] <- min(data) * 0.999
    }
    if (max(data) >= b_work[length(b_work)]) {
      b_work[length(b_work)] <- max(data) * 1.001
    }
    
    tryCatch({
      hist(data, breaks = b_work, plot = FALSE, right = FALSE, include.lowest = TRUE)$counts
    }, error = function(e) {
      rep(0, length(boundaries) - 1)
    })
  }
  
  d_vals <- numeric(reps)
  successful_reps <- 0
  
  for (i in seq_len(reps)) {
    # 1. Genera campione bootstrap
    N_tail_boot <- N - N_below
    
    if (N_tail_boot <= 0) {
      d_vals[i] <- NA_real_
      next
    }
    
    # Campioni sotto bmin (da distribuzione empirica)
    samp_below <- if (N_below > 0) {
      sample(x_pool_below, size = N_below, replace = TRUE)
    } else {
      numeric(0)
    }
    
    # Campioni sopra bmin (da power law teorica)
    samp_tail <- rpl_continuous(N_tail_boot, bmin_orig, alpha_orig)
    
    # Combina
    samp_all <- c(samp_below, samp_tail)
    
    # 2. Binna i dati bootstrap
    h_boot <- safe_binning(samp_all, boundaries)
    
    if (sum(h_boot) < N * 0.8) {  # Controllo perdita dati
      d_vals[i] <- NA_real_
      next
    }
    
    # 3. Estrai la porzione >= bmin 
    h_tail_boot <- h_boot[idx_bmin:length(h_boot)]
    
    # Assicura che h_tail_boot abbia la stessa lunghezza di h_tail
    if (length(h_tail_boot) != length(h_tail)) {
      # Tronca o estendi per avere la stessa lunghezza
      target_length <- length(h_tail)
      if (length(h_tail_boot) > target_length) {
        h_tail_boot <- h_tail_boot[1:target_length]
      } else if (length(h_tail_boot) < target_length) {
        # Estendi con zeri
        h_tail_boot <- c(h_tail_boot, rep(0, target_length - length(h_tail_boot)))
      }
    }
    
    if (sum(h_tail_boot) < 2) {
      d_vals[i] <- NA_real_
      next
    }
    
    n_tail_boot <- sum(h_tail_boot)
    
    # 4. Calcola statistica KS sui dati bootstrap, parametri originali (bmin_orig, alpha_orig)
    # CCDF empirica bootstrap
    ccdf_empirical_boot <- rev(cumsum(rev(h_tail_boot))) / n_tail_boot
    
    # CCDF teorica bootstrap 
    # Usa gli stessi bin_lowers_tail dei dati originali
    ccdf_theoretical_boot <- (bin_lowers_tail / bmin_orig)^(1 - alpha_orig)
    ccdf_theoretical_boot <- pmax(0, pmin(1, ccdf_theoretical_boot))
    
    # Verifica che abbiano la stessa lunghezza
    if (length(ccdf_empirical_boot) != length(ccdf_theoretical_boot)) {
      d_vals[i] <- NA_real_
      next
    }
    
    # Statistica KS bootstrap
    D_boot <- max(abs(ccdf_empirical_boot - ccdf_theoretical_boot))
    
    if (is.finite(D_boot)) {
      d_vals[i] <- D_boot
      successful_reps <- successful_reps + 1
    } else {
      d_vals[i] <- NA_real_
    }
    
    # Debug output
    if (!silent && (i <= 5 || i %% 200 == 0)) {
      cat(sprintf("Rep %d: D = %.6f (N_boot=%d)\n", i, D_boot, sum(h_boot)))
    }
  }
  
  # ===================================================================
  # CALCOLA P-VALUE
  # ===================================================================
  d_vals_clean <- d_vals[!is.na(d_vals) & is.finite(d_vals)]
  
  if (length(d_vals_clean) < 10) {
    warning("Troppo poche repliche bootstrap valide.")
    p_val <- NA_real_
  } else {
    # P-value = frazione di D_bootstrap >= D*
    p_val <- mean(d_vals_clean >= Dstar)
  }
  
  # ===================================================================
  # OUTPUT
  # ===================================================================
  if (!silent) {
    cat(sprintf("\nBootstrap completato:\n"))
    cat(sprintf("Repliche valide: %d/%d (%.1f%%)\n", 
                length(d_vals_clean), reps, 100*length(d_vals_clean)/reps))
    cat(sprintf("D* = %.6f\n", Dstar))
    
    if (length(d_vals_clean) > 0) {
      cat(sprintf("D bootstrap: media=%.6f, range=[%.6f, %.6f]\n", 
                  mean(d_vals_clean), min(d_vals_clean), max(d_vals_clean)))
      cat(sprintf("P-value = %.4f\n", p_val))
      
      cat(sprintf("D* vs media D_boot: %.6f vs %.6f\n", Dstar, mean(d_vals_clean)))
      cat(sprintf("Frazione D_boot >= D*: %.4f\n", mean(d_vals_clean >= Dstar)))
    }
  }
  
  return(list(
    p = p_val,
    d_vals = d_vals,
    Dstar = Dstar,
    n_valid_bootstrap = length(d_vals_clean),
    n_tail = n_tail
  ))
}



