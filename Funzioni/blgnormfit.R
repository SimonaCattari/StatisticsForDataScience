# La funzione stima i parametri di una distribuzione lognormale (mu e sigma) a partire da dati empirici
# Si considerano solo i bin in coda alla distribuzione
# Tramite una prima grid search, si valuta la log-likelihood per tutte le combinazioni di mu e sigma
# Se fine = TRUE, viene eseguita un'ulteriore grid search centrata sui valori migliori trovati in precedenza
# La funzione restituisce i parametri stimati e la log-likelihood massima


blgnormfit <- function(h, boundaries, murng = NULL, sigrng = NULL, bmin = NULL, fine = TRUE, fine_params = c(3, 1, 0.05)) {
  
  # Controlli iniziali
  if (!all(h == floor(h)) || any(h < 0)) stop("h must be a non-negative integer")
  if (length(boundaries) != (length(h) + 1)) stop("n boundaries must be n bins + 1")
  if (length(h) < 2) stop("at least 2 bins required")
  
  # Range di default 
  if (is.null(murng)) murng <- seq(-100, 50, 1)
  if (is.null(sigrng)) sigrng <- seq(1, 20, 0.1)
  
  # Gestione bmin
  if (is.null(bmin)) bmin <- boundaries[1]
  bmin_idx <- max(which(boundaries <= bmin))
  bmin <- boundaries[bmin_idx]
  
  # Seleziona dati >= bmin
  h2 <- h[bmin_idx:length(h)]
  l <- boundaries[bmin_idx:(length(boundaries) - 1)]
  u <- boundaries[(bmin_idx + 1):length(boundaries)]
  
  # Log-likelihood
  loglik_grid <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    
    # Probabilità teoriche per ogni bin
    p <- pnorm(log(u), mean = mu, sd = sigma) - pnorm(log(l), mean = mu, sd = sigma)
    
    # Normalizzazione per la coda >= bmin
    p <- p / (1 - pnorm(log(bmin), mean = mu, sd = sigma))
    
    if (any(p <= 0) || any(!is.finite(p))) return(-Inf)
    
    return(sum(h2 * log(p)))
  }
  
  # Grid search 
  ll_matrix <- outer(murng, sigrng, Vectorize(function(mu, sigma) loglik_grid(mu, sigma)))
  cat("Grid search: max loglik =", max(ll_matrix), "\n")
  
  # Trova l'indice del massimo
  max_idx <- which(ll_matrix == max(ll_matrix), arr.ind = TRUE)
  
  # 
  if (nrow(max_idx) > 0) {
    coarse_mu <- murng[max_idx[1, 1]]      
    coarse_sigma <- sigrng[max_idx[1, 2]]  
  } else {
    stop("No max found on grid search")
  }
  
  cat("Coarse estimates: mu =", coarse_mu, ", sigma =", coarse_sigma, "\n")
  
  if (fine) {
    fine_mu <- seq(coarse_mu - fine_params[1], coarse_mu + fine_params[1], by = fine_params[3])
    fine_sigma <- seq(max(0.01, coarse_sigma - fine_params[2]), coarse_sigma + fine_params[2], by = fine_params[3])
    # LL matrix per tutte le combinazioni mu-sigma
    ll_matrix_fine <- outer(fine_mu, fine_sigma, Vectorize(function(mu, sigma) loglik_grid(mu, sigma)))
    max_idx_fine <- which(ll_matrix_fine == max(ll_matrix_fine), arr.ind = TRUE)
    
    # Se necessario, aggiorna le stime 
    if (nrow(max_idx_fine) > 0) {
      mu_est <- fine_mu[max_idx_fine[1, 1]]
      sigma_est <- fine_sigma[max_idx_fine[1, 2]]
      L <- -max(ll_matrix_fine)
    } else {
    # Altrimenti utilizza i valori coarse
      mu_est <- coarse_mu
      sigma_est <- coarse_sigma
      L <- -max(ll_matrix)
    }
  } else { # Se fine = FALSE, restituisce le stime della grid search iniziale
    mu_est <- coarse_mu
    sigma_est <- coarse_sigma
    L <- -max(ll_matrix)
  }
  
  return(list(mu = mu_est, sigma = sigma_est, logLik = -L, bmin = bmin))
}

# ALTERNATIVA: ottimizzazione numerica invece di grid search
blgnormfit_optim <- function(h, boundaries, bmin = NULL) {
  if (!all(h == floor(h)) || any(h < 0)) stop("h must be a non-negative integer")
  if (length(boundaries) != (length(h) + 1)) stop("n boundaries must be n bins + 1")
  if (length(h) < 2) stop("at least 2 bins required")
  
  if (is.null(bmin)) bmin <- boundaries[1]
  bmin_idx <- max(which(boundaries <= bmin))
  bmin <- boundaries[bmin_idx]
  
  h2 <- h[bmin_idx:length(h)]
  l <- boundaries[bmin_idx:(length(boundaries) - 1)] # upper boundaries
  u <- boundaries[(bmin_idx + 1):length(boundaries)] # lower boundaries
  
  # Negative log-likelihood per optim()
  neg_loglik <- function(params) {
    mu <- params[1]
    sigma <- params[2]
    
    if (sigma <= 0) return(1e10) # Penalizzazione per sigma non valido
    
    # Probabilità che un valore cada all'interno di un intervallo
    p <- pnorm(log(u), mean = mu, sd = sigma) - pnorm(log(l), mean = mu, sd = sigma)
    # Normalizzazione rispetto a bmin
    p <- p / (1 - pnorm(log(bmin), mean = mu, sd = sigma))
    
    if (any(p <= 0) || any(!is.finite(p))) return(1e10)
    
    return(-sum(h2 * log(p)))
  }
  
  # Ottimizzazione numerica con starting values 
  # Stima iniziale: media dei log dei centri dei bin pesata per le frequenze
  bin_centers <- sqrt(l * u)
  mu_start <- sum(h2 * log(bin_centers)) / sum(h2)
  sigma_start <- 1.0
  
  opt <- optim(par = c(mu_start, sigma_start), fn = neg_loglik, method = "Nelder-Mead")
  
  return(list(mu = opt$par[1], sigma = opt$par[2], logLik = -opt$value, bmin = bmin))
}