# Power-Law Distributions in Binned Empirical Data: un'implementazione in R

Progetto per il corso di **Statistics for Data Science** (A.A. 2024/2025), basato sul paper *"Power-law distributions in binned empirical data"* di Virkar e Clauset (2014).

## Panoramica

Molti fenomeni naturali e sociali (terremoti, traffico sui siti web, ricchezza, dimensione delle città, dimensione degli incendi) vengono spesso descritti tramite distribuzioni power-law. Tuttavia, i dati disponibili sono spesso forniti in forma **binned** (istogrammi), il che rende inefficaci i classici metodi di stima e può portare a stime distorte, perdita di informazione e falsi positivi nel riconoscimento di una power-law.

Il progetto implementa interamente in **R** il framework statistico proposto da Virkar e Clauset, basato su:

- **Maximum Likelihood Estimation (MLE)** per la stima dei parametri
- **Statistica di Kolmogorov–Smirnov (KS)** per la selezione del limite inferiore `b_min`
- **Bootstrap semi-parametrico** per il test di bontà del fit (goodness-of-fit)
- **Test di Vuong** per il confronto tra modelli alternativi

L'obiettivo è fornire un toolkit affidabile per: stimare accuratamente i parametri di una power-law da dati a istogramma, testarne la plausibilità statistica ed evitare falsi positivi causati dagli artefatti introdotti dal binning.

## Metodologia

Il workflow implementato si articola in tre passaggi principali:

1. **Fit** — stima dei parametri ottimali (esponente di scala `α` e soglia `b_min`)
2. **Test** — valutazione della bontà del fit tramite p-value (modello accettato se `p ≥ 0.1`)
3. **Compare** — confronto del modello power-law con distribuzioni alternative (log-normale, esponenziale, esponenziale allungata, power-law con cutoff)

L'intero processo è stato validato sia su **dati sintetici** generati da distribuzioni note, sia su **dataset reali**, utilizzando diversi schemi di binning (lineare e logaritmico).

## Funzioni implementate

Sono state sviluppate 8 funzioni R principali:

1. bplfit(h, boundaries)  
Scopo: Stima i parametri della distribuzione power law per dati binned.  
Metodologia Statistica:
- Stima di α: Maximum Likelihood Estimation (MLE)
- Selezione di bmin: Minimizzazione statistica KS
- Ottimizzazione: Grid search su bmin, raffinamento numerico per α  
Output: alpha, bmin, logLik, D (statistica KS)  
Note: Implementazione dell'Algoritmo 1 del paper, gestione edge cases, controlli di convergenza.

2. getPDF(type, boundaries, bmin, ...)  
Scopo: Calcola la PDF per diverse distribuzioni su bins discreti.  
Distribuzioni supportate:
- Power law
- Exponential
- Log-normal
- Stretched exponential
- Power law con cutoff  
Metodologia: Integrazione delle PDF continue sui bin secondo la formula del paper (eq. 2.4).  
Normalizzazione garantita.

3. bplpva(h, boundaries, alpha, bmin, ntrials = 1000)  
Scopo: Test di goodness-of-fit per power law (bootstrap semi-parametrico).  
Metodologia:
- Generazione bootstrap: power law sopra bmin, resample sotto
- Statistica KS su ogni replica
- Calcolo p-value: P(D >= D*)  
Criterio di rigetto: p < 0.1  
Note: Re-stima bmin e α in ogni replica, uso inversione per generazione power law.

4. blrtest(h, boundaries, type1, params1, type2, params2)  
Scopo: Likelihood Ratio Test per modelli alternativi.  
Metodo di Vuong (1989):
- Statistica: R = somma log-likelihood ratio bin-wise
- Normalizzazione: R̃ = R / sqrt(n * var)
- P-value: H₀ = modelli equivalenti, H₁ = uno è superiore  
Interpretazione:
- R̃ > 0: favore a modello 1
- R̃ < 0: favore a modello 2
- |R̃| > 2 e p < 0.05: evidenza forte

**Funzioni di Fitting Alternative**

bexpnfit(h, boundaries, bmin)  
- Stima λ (esponenziale) via MLE  
- Ottimizzazione: Nelder-Mead + BFGS

blgnormfit(h, boundaries, bmin, murng, sigrng, fine = FALSE)  
- Stima μ e σ (log-normale)  
- Grid search + raffinamento fine  
- Controlli overflow/underflow

bstexpfit(h, boundaries, bmin)  
- Stima λ e β (stretched exp.)  
- Grid search bidimensionale + ottimizzazione numerica robusta

bplcutfit(h, boundaries, bmin)  
- Stima α e λ (power law con cutoff)  
- Ottimizzazione congiunta (multi-start, gradient-based)
  
Ogni funzione è stata validata tramite una batteria di test mirati (robustezza al rumore, effetto della dimensione campionaria, accuratezza su diversi schemi di binning, gestione di casi limite).

## Risultati principali

### Dati sintetici
Su 10.000 valori generati da una power-law nota (`α = 2.5`): stima `α = 2.512` (errore 0.4%), `p-value = 0.93`, con corretto rigetto dei modelli alternativi.

### Dati reali — Dimensione delle città USA
Dataset di 19.447 città statunitensi. Risultati confrontati con quelli del paper originale:

| Parametro | Nostra stima | Paper | Match |
|---|---|---|---|
| α | 2.376 | 2.38 | ✅ |
| b_min | 65.536 | 65.536 | ✅ |
| support | Moderate | Moderate | ✅ |

### Dati reali — Dimensione degli incendi (acri)
Dataset di 203.785 incendi su terreni federali USA (1986-1996):

| Parametro | Nostra stima | Paper | Match |
|---|---|---|---|
| α | 1.482 | 1.482 | ✅ |
| b_min | 2 | 2 | ✅ |
| support | None | None | ✅ |

In questo caso il modello power-law puro viene rigettato (`p-value = 0.000`); il modello che meglio si adatta ai dati è la **power-law con cutoff**.

## Conclusioni

- Il toolkit R sviluppato replica fedelmente i risultati del paper originale su entrambi i dataset reali analizzati.
- I dati binned introducono sfide specifiche (perdita di dettaglio, pattern nascosti, stime distorte, maggiore incertezza), che richiedono metodi statistici dedicati e spesso più dati per ottenere risultati affidabili.
- La scelta del metodo statistico ha un impatto significativo sui risultati ottenuti.
- Le power-law "pure" sono rare in natura: molti fenomeni si adattano meglio a modelli alternativi (es. power-law con cutoff), evidenziando l'importanza del confronto sistematico tra modelli.

## Dataset utilizzati

- **U.S. city sizes** — popolazione delle città statunitensi, censimento USA 2000
- **Wildfire sizes** — dimensione in acri degli incendi su terreni federali USA, 1986-1996

## Riferimenti

- Y. Virkar, A. Clauset. *Power-law distributions in binned empirical data*. The Annals of Applied Statistics, 2014.

## Struttura del repository

```
progetto13/
├── funzioni/
│   ├── bplfit.R          # Stima dei parametri della power-law (MLE + KS)
│   ├── getPDF.R          # PDF teorica su dati binned (5 modelli)
│   ├── bplpva.R          # Test di bontà del fit (bootstrap semi-parametrico)
│   ├── blrtest.R         # Confronto tra modelli (test di Vuong)
│   ├── bexpnfit.R        # Fit modello esponenziale
│   ├── blgnormfit.R      # Fit modello log-normale
│   ├── bstexpfit.R       # Fit modello esponenziale allungata
│   └── bplcutfit.R       # Fit modello power-law con cutoff
├── tests/                # Test di validazione per ciascuna funzione
    ├── bplfit_test.R
│   ├── getPDF_test.R         
│   ├── bplcutfit_test.R          
│   ├── bplpva_test.R          
│   ├── blrtest_test.R         
│   ├── bexpnfit_test.R        
│   ├── blgnormfit_test.R      
│   ├── bstexpfit_test.R
│   ├── city_test.R
│   ├── fire_test.R       
│   └── bplcutfit.R
├── 13.VirkarClauset.pdf # Paper di riferimento
├── Presentation13Cattari-Toscano-Trivelli.pptx # Presentazione lavoro
└── README.md
```

> Adatta questa sezione alla struttura reale del repository.

## Autori

- Cattari Simona
- Toscano Giuseppe Vincenzo Dylan
- Trivelli Matteo
