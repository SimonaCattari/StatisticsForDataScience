# Stima dei parametri lambda di una distribuzione esponenziale sui dati binned,
# usando mle a partire da un certo bmin.
# Calcola la probabilità teorica in ogni bin con formula chiusa.
# Permette sia grid search su lambda_range, sia ottimizzazione numerica.
# Restituisce il valore stimato di lambda e la log-likelihood massima.

bexpnfit <- function(h, boundaries, lambda_range = NULL, bmin = NULL) {
  
  # Controlli preliminari
  if (!all(h == floor(h)) || any(h < 0)) stop("h must be a non-negative integer")
  if (length(boundaries) != (length(h) + 1)) stop("n boundaries must be n bins + 1")
  if (length(h) < 2) stop("at least 2 bins required")
  
  # Gestione bmin
  if (is.null(bmin)) bmin <- boundaries[1]
  bmin_idx <- max(which(boundaries <= bmin))
  bmin <- boundaries[bmin_idx]
  
  # bin ≥ bmin
  h2 <- h[bmin_idx:length(h)]
  l <- boundaries[bmin_idx:(length(boundaries) - 1)]
  u <- boundaries[(bmin_idx + 1):length(boundaries)]
  
  if (length(h2) < 2 || sum(h2) == 0) {
    stop("Insufficient data above bmin")
  }
  
  # negative ll
  neg_loglik <- function(lambda) {
    if (lambda <= 0) return(1e10)
    
    # p_i per ciascun bin (formula chiusa -> probabilità condizionata)
    p <- (exp(-lambda * l) - exp(-lambda * u)) / exp(-lambda * bmin)
    
    if (any(p <= 0) || any(!is.finite(p))) return(1e10)
    
    return(-sum(h2 * log(p)))
  }
  
  # Stima di lambda
  if (!is.null(lambda_range)) {
    nll_values <- sapply(lambda_range, neg_loglik)
    min_idx <- which.min(nll_values)
    lambda_est <- lambda_range[min_idx]
    L <- -nll_values[min_idx]
  } else {
    opt <- optim(par = 0.01, fn = neg_loglik, method = "Nelder-Mead", lower = 1e-6, upper = 5)
    lambda_est <- opt$par
    L <- -opt$value
  }
  
  return(list(lambda = lambda_est, loglik = L))
}
