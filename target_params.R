### --- Vaccination serology target parameters --- ###

# Use low concentration S in serology
low_conc_s <- F

# Cluster colours
# For vaccination/infection "clusters" of the questionnaire
clust_cols <- RColorBrewer::brewer.pal(9, "Purples")[seq(2, 8, 2)]
# For seroclusters
clust_cols2 <- c("#A9B4C2", "#7989A4", "#60718E", "#445D7B")

# Unified colours for the Abs, similar colours for related targets
# NOTE: using R pipe as the magrittr pipe would require loading the package outside of tar_option_set
ab_cols <- c(RColorBrewer::brewer.pal(9, "BuGn")[4:6] |>
               setNames(c("Anti-RBD", "Anti-S1", "Anti-S1S2")),
             RColorBrewer::brewer.pal(9, "Oranges")[c(3, 5)] |>
               setNames(c("Anti-Na", "Anti-Nc")),
             RColorBrewer::brewer.pal(9, "BuGn")[4] |>
               setNames("Anti-S"),
             RColorBrewer::brewer.pal(9, "Oranges")[4] |>
               setNames("Anti-N"))

# Linearising and scaling data before protein clustering
prot_lin <- T
prot_scale <- T