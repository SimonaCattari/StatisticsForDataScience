# bstexpfit: stima i parametri (lambda, beta) di una distribuzione stretched exponential 
# su dati binned tramite mle.
# Supporta sia una grid search su un range specificato, sia ottimizzazione numerica (Nelder-Mead).
# Il modello è troncato inferiormente a bmin e normalizzato di conseguenza.
# Restituisce le stime dei parametri e la log-likelihood massima.


bstexpfit <- function(h, boundaries, range = NULL, bmin = NULL) {
  
  if (!all(h == floor(h)) || any(h < 0)) stop("h must be a non-negative integer")
  if (length(boundaries) != (length(h) + 1)) stop("n boundaries must be n bins + 1")
  if (length(h) < 2) stop("at least 2 bins required")
  
  if (is.null(bmin)) {
    bmin <- boundaries[1]
  }
  
  # Trova il primo indice per cui 
  ind <- max(which(boundaries <= bmin))
  bmin <- boundaries[ind]
  
  h2 <- h[ind:length(h)]
  boundaries2 <- boundaries[ind:(ind + length(h2))]
  
  l <- boundaries2[-length(boundaries2)]
  u <- boundaries2[-1]
  
  compute_loglik <- function(lambda, beta) {
    if (lambda <= 0 || beta <= 0) return(-Inf)
    
    # Probabilità teoriche per ogni bin [l, u)
    # Per stretched exponential: f(x) = x^(β-1) * exp(-λx^β)
    probs <- exp(-lambda * l^beta) - exp(-lambda * u^beta)
    
    # Costante di normalizzazione: P(X >= bmin)
    norm_const <- exp(-lambda * bmin^beta)
    
    # Normalizza le probabilità
    probs <- probs / norm_const
    
    # Evita log(0)
    probs[probs <= 0] <- 1e-10
    
    # Log-likelihood 
    logL <- sum(h2 * log(probs))
    return(logL)
  }
  
  if (!is.null(range)) {
    # Grid search 
    grid <- expand.grid(lambda = range$lambda, beta = range$beta)
    logliks <- apply(grid, 1, function(p) compute_loglik(p[1], p[2]))
    best_index <- which.max(logliks)
    
    return(list(
      lambda = grid$lambda[best_index],
      beta = grid$beta[best_index],
      loglik = logliks[best_index]
    ))
    
  } else {
    # Ottimizzazione numerica
    negloglik <- function(params) {
      lambda <- params[1]; beta <- params[2]
      return(-compute_loglik(lambda, beta))
    }
    
    # Usa Nelder-Mead per problemi non-derivabili
    opt <- optim(c(0.1, 0.5), negloglik, method = "Nelder-Mead", 
                 control = list(maxit = 10000))
    
    if (opt$convergence != 0) {
      warning("Optimization may not have converged")
    }
    
    return(list(
      lambda = opt$par[1], 
      beta = opt$par[2], 
      loglik = -opt$value
    ))
  }
}