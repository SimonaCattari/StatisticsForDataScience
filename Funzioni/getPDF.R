# ===================================================================
# getPDF.R: Prob che il valore cada all'interno di un bin.
#           Calcola una distribuzione di probabilità (PDF) per dati
#           binned per i seguenti modelli:
#             - pl (power law);
#             - expn (esponenziale)
#             - lgnorm (lognormale)
#             - stexp (stretched exponential)
#             - plcut (power law con cutoff esponenziale)
# Restituisce un vettore di probabilità normalizzate (somma = 1), una
# per ciascun bin compreso tra il limiti superiori e inferiori di 
# interesse
# ===================================================================

getPDF <- function(boundaries, type, bmin = NULL, ...) {
  args <- list(...)
  
  # Almeno 3 estremi di bin per costruitr almeno 2 intervalli
  if (length(boundaries) < 3) {
    cat("Error: at least 3 bin boundaries are required\n")
    return(NULL)
  }
  
  # Validazione bmin: deve essere un numero finito singolo, maggiore del limite minimo
  if (is.null(bmin) || !is.numeric(bmin) || length(bmin) != 1 || !is.finite(bmin) || bmin < boundaries[1]) {
    bmin <- boundaries[1]
    warning("Invalid or missing bmin, defaulting to boundaries[1].")
  }
  
  # Indice del primo bin valido
  ind <- max(which(boundaries <= bmin))
  if (length(ind) == 0 || is.na(ind)) {
    ind <- 1 # default
  }
  
  # Definisco i limiti inferiori e superiori dei bin sopra bmin
  l_relevant <- boundaries[ind:(length(boundaries) - 1)]
  u_relevant <- boundaries[(ind + 1):length(boundaries)]
  
  # Inizializza probabilità non normalizzate e costante di normalizzazione
  prob_unnorm_segment <- NULL
  normalization_const <- 1
  
  # Calcolo PDF in base al tipo richiesto per ciascun bin
  prden <- switch(type,
                  
                  # Power law
                  pl = {
                    alpha <- args[[1]]
                    if (is.null(alpha) || !is.numeric(alpha) || !is.finite(alpha) || alpha <= 1) {
                      cat("Error: alpha must be > 1 for power-law\n")
                      return(NULL)
                    }
                    
                    normalization_const <- bmin^(1-alpha)
                    prob_unnorm_segment <- l_relevant^(1-alpha) - u_relevant^(1-alpha)
                    prob_unnorm_segment[prob_unnorm_segment <= 0] <- .Machine$double.xmin
                    prob_unnorm_segment
                  },
                  
                  # Exponential
                  expn = {
                    lambda <- args[[1]]
                    if (is.null(lambda) || !is.numeric(lambda) || !is.finite(lambda) || lambda <= 0) {
                      cat("Error: lambda must be > 0 for exponential\n")
                      return(NULL)
                    }
                    
                    normalization_const <- exp(-lambda * bmin)
                    
                    # Caso di underflow
                    prob_unnorm_segment <- tryCatch({
                      exp_l <- exp(-lambda * l_relevant)
                      exp_u <- exp(-lambda * u_relevant)
                      
                      # Se i valori sono troppo piccoli, usa log-differenze per evitare underflow
                      if (any(exp_l < .Machine$double.xmin) || any(exp_u < .Machine$double.xmin)) {
                        log_diff <- pmax(-lambda * l_relevant, log(.Machine$double.xmin)) - 
                          pmax(-lambda * u_relevant, log(.Machine$double.xmin))
                        exp(log_diff)
                      } else {
                        exp_l - exp_u
                      }
                    }, error = function(e) {
                      warning("Exponential PDF: overflow occurred, using uniform distribution instead")
                      rep(1/length(l_relevant), length(l_relevant))
                    })
                    
                    prob_unnorm_segment[prob_unnorm_segment <= 0] <- .Machine$double.xmin
                    prob_unnorm_segment
                  },
                  
                  # Log-normal
                  lgnorm = {
                    mu <- args[[1]]
                    sigma <- args[[2]]
                    if (is.null(mu) || is.null(sigma) || !is.numeric(mu) || !is.numeric(sigma) ||
                        !is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
                      cat("Error: mu and sigma must be finite and sigma > 0 for log-normal\n")
                      return(NULL)
                    }
                    
                    # Gestione casi estremi per mu e sigma
                    if (abs(mu) > 50 || sigma > 10) {
                      warning("Extreme log-normal parameters; results may be unstable")
                    }
                    
                    normalization_const <- tryCatch({
                      1 - pnorm(log(bmin), mean = mu, sd = sigma)
                    }, error = function(e) {
                      .Machine$double.xmin
                    })
                    
                    # Se la normalizzazione è troppo piccola, si usa una distribuzione uniforme
                    if (normalization_const <= .Machine$double.xmin) {
                      warning("Log-normal normalization too small, using uniform distribution")
                      return(rep(1/length(l_relevant), length(l_relevant)))
                    }
                    
                    # Massa di probabilità in ciascun bin (approssimazione tramite CDF lognormale)
                    prob_unnorm_segment <- tryCatch({
                      pnorm(log(u_relevant), mean = mu, sd = sigma) - 
                        pnorm(log(l_relevant), mean = mu, sd = sigma)
                    }, error = function(e) {
                      warning("Log-normal PDF: error during computation, using uniform distribution")
                      rep(1/length(l_relevant), length(l_relevant))
                    })
                    
                    prob_unnorm_segment[prob_unnorm_segment <= 0] <- .Machine$double.xmin
                    prob_unnorm_segment
                  },
                  
                  # Stretched exponential
                  stexp = {
                    lambda <- args[[1]]
                    beta <- args[[2]]
                    if (is.null(lambda) || is.null(beta) || !is.numeric(lambda) || !is.numeric(beta) ||
                        !is.finite(lambda) || !is.finite(beta) || lambda <= 0 || beta <= 0) {
                      cat("Error: lambda > 0 and beta > 0 required for stretched exponential\n")
                      return(NULL)
                    }
                    
                    normalization_const <- exp(-lambda * bmin^beta)
                    
                    # Gestione overflow
                    prob_unnorm_segment <- tryCatch({
                      term_l <- -lambda * l_relevant^beta
                      term_u <- -lambda * u_relevant^beta
                      term_l[term_l < -700] <- -700 
                      term_u[term_u < -700] <- -700
                      
                      exp(term_l) - exp(term_u)
                    }, error = function(e) {
                      warning("Stretched exponential: overflow gestito")
                      rep(1/length(l_relevant), length(l_relevant))
                    })
                    
                    prob_unnorm_segment[prob_unnorm_segment <= 0] <- .Machine$double.xmin
                    prob_unnorm_segment
                  },
                  
                  # Power law con cutoff esponenziale
                  plcut = {
                    alpha <- args[[1]]
                    lambda <- args[[2]]
                    if (is.null(alpha) || is.null(lambda) || !is.numeric(alpha) || !is.numeric(lambda) ||
                        !is.finite(alpha) || !is.finite(lambda) || alpha <= 1 || lambda <= 0) {
                      cat("Error: alpha > 1 and lambda > 0 required for power-law with cutoff\n")
                      return(NULL)
                    }
                    
                    # Approssimazione più stabile per PL+cutoff
                    prob_unnorm_segment <- tryCatch({
                      
                      pdf_l <- l_relevant^(-alpha) * exp(-lambda * l_relevant)
                      pdf_u <- u_relevant^(-alpha) * exp(-lambda * u_relevant)
                      
                      # Controlla overflow
                      pdf_l[!is.finite(pdf_l)] <- .Machine$double.xmin
                      pdf_u[!is.finite(pdf_u)] <- .Machine$double.xmin
                      
                      (pdf_l + pdf_u) / 2 * (u_relevant - l_relevant) # regola del trapezio perchè non si può usare l'integrale esatta
                    }, error = function(e) {
                      warning("Power-law with cutoff: computation failed, using uniform distribution.")
                      rep(1/length(l_relevant), length(l_relevant))
                    })
                    
                    prob_unnorm_segment[prob_unnorm_segment <= 0] <- .Machine$double.xmin
                    prob_unnorm_segment
                  },
                  
                  {
                    cat("Error: unknown distribution type:", type, "\n")
                    return(NULL)
                  }
  )
  
  # Se c'è errore, interrompe
  if (is.null(prob_unnorm_segment)) return(NULL)
  
  # Normalizzazione iniziale: divisione delle probabilità non normalizzate per la costante
  prden_raw <- tryCatch({
    prob_unnorm_segment / normalization_const
  }, error = function(e) {
    warning("Error during initial normalization.")
    prob_unnorm_segment
  })
  
  # Controllo che la densità sia valida
  if (any(!is.finite(prden_raw)) || any(is.na(prden_raw)) || sum(prden_raw) <= 0) {
    warning("Invalid probabilities after normalization, using uniform distribution.")
    return(rep(1/length(l_relevant), length(l_relevant))) # stessa probabilità a tutti i bin
  }
  
  # Normalizzazione 
  prden_final <- prden_raw / sum(prden_raw)
  
  # Check finale
  if (any(!is.finite(prden_final)) || any(prden_final < 0)) {
    warning("Final probability check failed, using uniform distribution.")
    return(rep(1/length(l_relevant), length(l_relevant)))
  }
  
  return(as.numeric(prden_final))
}