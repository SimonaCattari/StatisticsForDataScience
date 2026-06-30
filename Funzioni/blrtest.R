# =============================================================================
# blrtest(): Esegue un test del rapporto di verosimiglianza (Likelihood Ratio Test)
#            tra due modelli alternativi p1 e p2 per dati binned.
#
# USO:
#   Confronta due distribuzioni di probabilità teoriche `p1` e `p2` (vettori di 
#   probabilità normalizzate per ciascun bin), rispetto a una distribuzione empirica
#   binned `h` con confini `boundaries`, considerando solo la parte di dati >= bmin.
#
# ARGOMENTI:
#   - p1, p2:    vettori di probabilità teorica (una per bin), normalizzati
#   - h:         vettore di conteggi per ciascun bin
#   - boundaries: estremi dei bin (length = length(h) + 1)
#   - bmin:      soglia minima (es. valore stimato di bmin)
#   - isNested:  se TRUE, assume che i modelli siano nidificati (usa test semplificato)
#
# OUTPUT:
#   - normR:     valore normalizzato della statistica del rapporto di verosimiglianza
#   - p:         p-value del test (probabilità che R ≥ osservato sotto H0)
#   - R:         valore grezzo del log-likelihood ratio
#   - n:         numero di osservazioni nella coda (≥ bmin)
# =============================================================================


blrtest <- function(p1, p2, h, boundaries, bmin, isNested = FALSE) {
  
  # h deve essere un vettore di conteggi interi non negativi
  if (!all(h == floor(h))) stop("Vector h should be an integer vector.")
  if (any(h < 0)) stop("Vector h should be non-negative.")
  
  # boundaries deve avere una lunghezza pari a length(h) + 1
  if (length(boundaries) != (length(h) + 1)) stop("Incorrect number of elements in either boundaries or h.")
  if (length(h) < 2) stop("Need at least 2 bins to make this work.")
  
  # Converto tutti gli input in vettori numerici
  h <- as.vector(h)
  boundaries <- as.vector(boundaries)
  p1 <- pmax(as.vector(p1), 1e-10)  # Evito log(0)
  p2 <- pmax(as.vector(p2), 1e-10)
  
  # Trova l'indice del primo bin <= bmin
  if (all(boundaries > bmin)) stop("bmin is less than all boundary values.")
  ind <- max(which(boundaries <= bmin))
  
  # =================================
  # ESTRAI I DATI NELLA CODA (≥ bmin)
  # =================================
  
  h2 <- h[ind:length(h)] # conteggi dei bin nella coda
  boundaries2 <- boundaries[ind:(ind + length(h2))] # confini corrispondenti
  
  l <- boundaries2[-length(boundaries2)]
  u <- boundaries2[-1]
  temp <- (l + u) / 2 # centroidi di ciascun bin (media tra l e u)
  
  if (length(temp) != length(h2)) stop("Mismatch between temp and h2 lengths.")
  
  
  temp2 <- unlist(mapply(rep, temp, h2)) # associa ogni osservazione simulata a un valore continuo (necessario per calcolare la likelihood)
  
  n <- sum(h2) # osservazioni totali nella coda
  N <- length(temp2)
  set.seed(42)
  temp2 <- sample(temp2, N) # Permuta casualmente l'ordine
  
  whichbin <- findInterval(temp2, boundaries2, rightmost.closed = TRUE) # assrgna ciascuna osservazione al bin corrispondente
  if (any(whichbin < 1 | whichbin > length(p1))) stop("Bin index out of bounds in p1/p2.")
  
  l1 <- log(p1[whichbin])
  l2 <- log(p2[whichbin])
  
  l1_bar <- mean(l1)
  l2_bar <- mean(l2)
  R <- sum(l1 - l2) # LLR
  
  # se un modello non è incluso nell'altro
  if (!isNested) {
    temp <- (l1 - l2 - (l1_bar - l2_bar))^2 # Varianza campionaria della differenza normalizzata
    sigmaR <- sqrt(mean(temp, na.rm = TRUE)) # Stima della sd di R
    if (sigmaR == 0) {
      warning("Standard deviation sigmaR is zero; returning NA.")
      normR <- NA
      p <- NA
    } else {
      normR <- R / (sqrt(n) * sigmaR) # statistica normalizzata
      p <- erfc(abs(R) / (sqrt(2 * n) * sigmaR)) # p-value
    }
  } else {
    p <- erfc(sqrt(abs(R)) / sqrt(2)) # p-value viene stimato direttamente
    normR <- R
  }
  
  return(list(normR = normR, p = p, R = R, n = n))
}

erfc <- function(x) 2 * pnorm(-x * sqrt(2)) # Trasforma R in un p-value
