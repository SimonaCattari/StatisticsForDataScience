# bplcutfit: stima i parametri (alpha, lambda) di una distribuzione power-law con cutoff esponenziale
# a partire da dati binned, usando mle.
# Supporta due modalità: grid search (se vengono forniti range di alpha e lambda) oppure ottimizzazione numerica.
# Il modello utilizza un troncamento inferiore in corrispondenza di bmin.
# Restituisce le stime dei parametri e la log-likelihood massima.


bplcutfit <- function(h, boundaries, range_alpha = NULL, range_lambda = NULL, bmin = NULL) {
  
  if (!all(h == floor(h)) || any(h < 0)) stop("h must be a non-negative integer")
  if (length(boundaries) != (length(h) + 1)) stop("n boundaries must be n bins + 1")
  if (length(h) < 2) stop("at least 2 bins required")
  
  
  if (is.null(bmin)) bmin <- boundaries[1]
  
  # Prendo il lower bound del bin in cui cade bmin
  bmin_idx <- max(which(boundaries <= bmin))
  bmin <- boundaries[bmin_idx]
  
  # Seleziono i dati sopra bmin
  h2 <- h[bmin_idx:length(h)]
  l <- boundaries[bmin_idx:(length(boundaries) - 1)] # lower bounds
  u <- boundaries[(bmin_idx + 1):length(boundaries)] # upper bounds
  
  if (length(h2) < 2 || sum(h2) == 0) {
    stop("Insufficient data above bmin")
  }
  
  # log-likelihood function
  loglik <- function(par){
    alpha <- par[1]
    lambda <- par[2]
    
    # Se i parametri sono fuori dal dominio valido, restituisce una likelihood molto bassa
    if (alpha <= 1 || lambda <= 0) return(-Inf)
    
    # Costante di normalizzazione 
    norm_const <- tryCatch(
      integrate(function(x) x^(-alpha) * exp(-lambda * x), lower = bmin, upper = Inf, rel.tol = 1e-8)$value,
      error = function(e) return(NA)
    )
    
    if (!is.finite(norm_const) || norm_const <= 0) return(-Inf)
    
    # Probabilità teorica di un valore di cadere in ciascun bin [l, u)
    p <- mapply(function(li, ui) {
      tryCatch(
        integrate(function(x) x^(-alpha) * exp(-lambda * x), lower = li, upper = ui, rel.tol = 1e-8)$value,
        error = function(e) return(0)
      )
    }, l, u)
    
    # Normalizzazione delle p
    p <- p / norm_const
    
    if (any(p <= 0) || any(!is.finite(p))) return(-Inf)
    
    # log-likelihood
    return(sum(h2 * log(p)))
  }
  
  # Grid search se sono forniti range di alpha e lambda
  if (!is.null(range_alpha) && !is.null(range_lambda)) {
    
    grid <- expand.grid(alpha = range_alpha, lambda = range_lambda)
    
    # Crea un vettore di log-likelihoods per ciascuna combinazione lambda-apha
    ll_values <- apply(grid, 1, function(par){
      loglik(c(par[1], par[2]))
    })
    
    # log-likelihood massima
    max_idx <- which.max(ll_values) 
    
    alpha_est <- grid$alpha[max_idx]
    lambda_est <- grid$lambda[max_idx]
    L <- ll_values[max_idx]
    
  } else {
    
    # Ottimizzazione numerica 
    neg_loglik <- function(par) -loglik(par)  # Converte in problema di minimizzazione 
    opt <- optim(par = c(2, 0.01), fn = neg_loglik, method = "Nelder-Mead")
    
    alpha_est <- opt$par[1]
    lambda_est <- opt$par[2]
    L <- -opt$value  # Riconverte la log-likelihood in positiva
  }
  
  return(list(alpha = alpha_est, lambda = lambda_est, loglik = L))
}