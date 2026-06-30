if (!exists("results_global")) {
  message("results_global non trovato: eseguo test_bplfit_complete() …")
  source("Test/bplfit_test.R")
  test_bplfit_complete() 
}

if (!requireNamespace("ggplot2", quietly = TRUE))   stop("Pacchetto ggplot2 mancante")
if (!requireNamespace("ggrepel",  quietly = TRUE))   stop("Pacchetto ggrepel mancante")

library(dplyr) # manipolazione dati
library(stringr) # stringhe
library(ggplot2)
library(ggrepel)

results_df <- results_global %>%
  mutate(family = case_when(
    str_detect(scenario, "Noisy")         ~ "Dati rumorosi",
    str_detect(scenario, "^N-")           ~ "Varia N",
    str_detect(scenario, "^Alpha-")       ~ "Varia α",
    str_detect(scenario, "Binning")       ~ "Schemi di binning",
    TRUE                                   ~ "Altro"
  ))

families_to_plot <- c("Dati rumorosi", "Varia N", "Varia α", "Schemi di binning")

# ---------------------------------------------------------------
# FUNZIONE DI PLOT
# ---------------------------------------------------------------
plot_family <- function(df, fam) {
  ggplot(df, aes(alpha_true, alpha_hat, label = scenario, size = bmin_error)) +
    geom_point(color = "#1F77B4") +
    geom_text_repel(size = 3, max.overlaps = Inf) +
    scale_size_continuous(name = "|Errore su b̂min|") +
    labs(title = paste0("α vero vs α stimato – ", fam),
         x = expression(alpha[true]),
         y = expression(hat(alpha))) +
    theme_minimal(base_size = 13)
}

# ---------------------------------------------------------------
# CREAZIONE E VISUALIZZAZIONE PLOT
# ---------------------------------------------------------------
plots <- lapply(families_to_plot, function(fam) {
  df_sub <- filter(results_df, family == fam, is.finite(bmin_error))
  if (nrow(df_sub) == 0) return(NULL)
  plot_family(df_sub, fam)
})

names(plots) <- families_to_plot

for (p in plots) if (!is.null(p)) print(p)
