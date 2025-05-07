# Proteomics target functions

# Correlate Olink and Alamar proteomics for overlapping binders
run_prot_corr <- function(dat_o, dat_a, outliers_o, outliers_a) {
  # Use "raw" data where proteins are not filtered out yet
  
  # Merge replicates
  dat_a <- combine_repl(dat_a)
  
  # Overlap in proteins and samples
  prot_overlap <- intersect(
    dat_a$binfo %>% select(assay, uniprot),
    dat_o$binfo %>% select(assay, uniprot)
  )
  
  # Exclude outliers as they affect the linear relation drastically
  excl_o <- outliers_o$data %>% filter(!is.na(outlier)) %>% pull(unique_sample_name)
  excl_a <- outliers_a$data %>% filter(!is.na(outlier)) %>% pull(unique_sample_name)
  samp_overlap <- intersect(
    dat_a$npm %>% filter(!unique_sample_name %in% excl_a) %>% pull(unique_sample_name),
    dat_o$npm %>% filter(!unique_sample_name %in% excl_o) %>% pull(unique_sample_name)
  )
  
  # Venn diagram
  venn_dat <- data.frame("protein" = unique(c(dat_a$binfo %>% pull(uniprot),
                                              dat_o$binfo %>% pull(uniprot)))) %>%
    mutate("Olink" = protein %in% (dat_o$binfo %>% pull(uniprot)),
           "Alamar" = protein %in% (dat_a$binfo %>% pull(uniprot)))
  gvenn <- ggvenn(venn_dat, columns = c("Olink", "Alamar"), text_size = 4)
  
  # Data for computing correlations
  comp_dat <- rbind(
    dat_a$npm %>% select(unique_sample_name, !!prot_overlap$assay) %>%
      filter(unique_sample_name %in% samp_overlap) %>%
      pivot_longer(cols = -unique_sample_name) %>% mutate(dataset = "alamar"),
    dat_o$npm %>% select(unique_sample_name, !!prot_overlap$assay) %>%
      filter(unique_sample_name %in% samp_overlap) %>%
      pivot_longer(cols = -unique_sample_name) %>% mutate(dataset = "olink")
  ) %>%
    pivot_wider(id_cols = c(unique_sample_name, name), names_from = "dataset", values_from = "value")
  
  comp_cor <- comp_dat %>%
    summarise(corr = cor(alamar, olink, method = "spearman"), .by = name)
  
  # Correlation with missing frequencies
  mf <- comp_cor %>%
    # Add missing frequencies (proportion of samples below LOD) to colour bars by datasets where protein is detected
    left_join(dat_a$binfo %>% select(assay, overall), by = c("name" = "assay")) %>%
    mutate(overall = (100 - overall) / 100) %>%
    rename(mf_a = overall) %>%
    left_join(dat_o$binfo %>% select(assay, missingfreq), by = c("name" = "assay")) %>%
    rename(mf_o = missingfreq) %>%
    mutate(det_which = case_when(mf_a > 0.5 & mf_o > 0.5 ~ "Not detected",
                                 mf_a < 0.5 & mf_o > 0.5 ~ "Detected Alamar",
                                 mf_a > 0.5 & mf_o < 0.5 ~ "Detected Olink",
                                 mf_a < 0.5 & mf_o < 0.5 ~ "Detected both")) %>%
    
    mutate(det_which = factor(det_which, levels = c("Detected both", "Detected Alamar", "Detected Olink", "Not detected"))) %>%
    # Bin data
    mutate(bin = cut(corr, breaks = seq(-1, 1, by = 0.2)))
  
  cor_barplt <- ggplot(mf %>% count(det_which, bin)) +
    geom_bar(aes(x = bin, y = n, fill = det_which), stat = "identity", show.legend = ifelse(length(unique(mf$det_which)) > 1, T, F)) +
    geom_text(data = mf %>% count(bin), mapping = aes(x = bin, y = n + 2, label = n), size = 6) +
    labs(x = "Correlation", y = "Count", fill = "Detection") +
    theme_classic(14)
  
  cor_scatterplt <- ggplot(mf, aes(x = 1 - mf_o, y = 1 - mf_a, colour = corr)) +
    geom_point(size = 3, alpha = 0.75) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey") +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey") +
    scale_colour_viridis_c() +
    labs(x = "Detectability PEA", y = "Detectability NULISA", colour = "Correlation") +
    theme_classic(14)
  
  # Per-protein scatter plots
  scatter_name <- make_fn("protein_correlation", "corr_scatter.pdf")
  
  pdf(scatter_name, onefile = T, height = 8, width = 8)
  par(mfrow = c(3, 3))
  
  for (prot in prot_overlap$assay) {
    plt_title <- paste0(prot, ", rho: ", comp_cor %>% filter(name == prot) %>% pull(corr) %>% signif(3))
    
    plot(
      dat_o$npm[match(samp_overlap, dat_o$npm$unique_sample_name), prot, drop = T] %>% scale(),
      dat_a$npm[match(samp_overlap, dat_a$npm$unique_sample_name), prot, drop = T] %>% scale(),
      xlab = "Olink", ylab = "Alamar", main = plt_title
    )
    abline(0, 1)
  }
  
  dev.off()
  
  return(list("n_prot" = nrow(prot_overlap), "n_samp" = length(samp_overlap),
              "corr" = mf, "barplt" = cor_barplt, "scatterplt" = cor_scatterplt,
              "prot_scatter_file" = scatter_name, "venn" = gvenn))
}

# Make plot showing number of proteins that are detected
make_prot_detectability_plot <- function(olink_dat, alamar_dat) {
  # Get detected percentage for all proteins in both data sets
  det_percent <- rbind(
    olink_dat$binfo %>% mutate(det_perc = (1 - missingfreq) * 100) %>% select(assay, det_perc),
    alamar_dat$binfo %>% select(assay, det_perc = overall)
  )
  
  plt_out <- ggplot(det_percent, aes(x = det_perc)) +
    geom_histogram(bins = 20, colour = "black", fill = "grey", breaks = seq(0, 100, 5)) +
    labs(x = "Protein detectability (% of samples)", y = "Protein count (N)") +
    theme_classic(16)
  
  return(list("data" = det_percent, "plot" = plt_out))
}

# Make some QC plots (sample outlier plots made in data processing)
# Protein CV plot looking at pooled samples
make_prot_cv_plot <- function(olink_raw, alamar_raw, prot_dat) {
  # Olink CV
  cv_olink <- olink_raw$npm_pool %>%
    pivot_longer(cols = -unique_sample_name) %>%
    summarise(cv = cv(2**value), .by = name) %>%
    left_join(olink_raw$binfo %>% select(assay, det_percent = missingfreq),
              by = c("name" = "assay")) %>%
    mutate(det_percent = (1 - det_percent) * 100) %>%
    mutate(method = "olink") %>%
    left_join(prot_dat$binfo %>% select(assay, unique_protein_name, method),
              by = c("name" = "assay", "method"))
  
  # Alamar CV
  alamar_pools <- alamar_raw$sinfo %>%
    filter(str_detect(unique_sample_name, "P0\\d_P\\d")) %>%
    pull(unique_sample_name)
  cv_alamar <- alamar_raw$npm %>%
    filter(unique_sample_name %in% alamar_pools) %>%
    pivot_longer(cols = -unique_sample_name) %>%
    summarise(cv = cv(value), .by = name) %>%
    left_join(alamar_raw$binfo %>% select(assay, det_percent = overall),
              by = c("name" = "assay")) %>%
    mutate(method = "alamar") %>%
    left_join(prot_dat$binfo %>% select(assay, unique_protein_name, method),
              by = c("name" = "assay", "method"))
  
  cvs <- rbind(cv_olink, cv_alamar)
  plt_out <- ggplot(cvs, aes(x = cv)) +
    geom_histogram(aes(fill = det_percent > 50), binwidth = 5, colour = "black",
                   breaks = seq(0, max(cvs$cv), 5)) +
    labs(x = "Coefficient of variation (%)", y = "Protein count (N)",
         fill = "> 50% detectability") +
    scale_fill_manual(values = c("TRUE" = "grey85", "FALSE" = "sienna2")) +
    theme_classic(16)
  
  return(list("data" = cvs, "plot" = plt_out))
}

# Sample PCA plot of PCs 1 and 2, coloured by questionnaire vaccination and infection status
make_prot_pca_plot <- function(prot_dat, clust_cols) {
  # Input data with annotation
  plt_dat <- prot_dat$npm %>%
    left_join(prot_dat$sinfo %>% select(unique_sample_name, vacc_inf),
              by = "unique_sample_name") %>%
    filter(!is.na(vacc_inf)) %>%
    # Set order of levels
    mutate(vacc_inf = factor(vacc_inf, levels = c("V- I-", "V+ I-", "V- I+", "V+ I+")))
  
  colrs <- clust_cols %>% setNames(levels(plt_dat$vacc_inf))
  
  plt_out <- plot_pca(plt_dat %>% select(-unique_sample_name, -vacc_inf), pcs = 1:2,
                      scale = T, center = T,
                      col_vec = pull(plt_dat, vacc_inf),
                      col_scale = scale_colour_manual(values = colrs)) +
    labs(colour = "Self-reported\nimmune status")
  
  return(plt_out)
}

# Plot principal components 1 to 4 in a paired plot, colouring by proteotype
make_prot_pca1to4 <- function(prot_dat, prot_clust) {
  
  plt_dat <- prot_dat$npm %>%
    left_join(prot_clust$clust_dat, by = "unique_sample_name") %>%
    mutate(cluster = as.character(cluster))
  
  colrs <- c("grey", brewer.pal(length(unique(plt_dat$cluster)) - 1, "Dark2"))
  
  plt_out <- plot_pca(plt_dat %>% select(-unique_sample_name, -cluster), pcs = 1:4,
                      col_vec = pull(plt_dat, cluster), scale = T, center = T,
                      col_scale = scale_colour_manual(values = colrs),
                      fill_scale = scale_fill_manual(values = colrs))
  return(plt_out)
}

# Protein UMAP plot
make_prot_umap_plot <- function(prot_in, rseed = 123) {
  
  if (is.numeric(rseed)) {set.seed(rseed)}
  
  # Run PCA
  ump <- prot_in$npm %>%
    column_to_rownames("unique_sample_name") %>%
    umap(n_neighbors = 10,
         a = 1.8956,
         b = 0.8006)
  
  # Make data for plotting
  plt_dat <- ump %>%
    as.data.frame() %>%
    select(V1, V2) %>%
    rownames_to_column("unique_sample_name") %>%
    left_join(prot_in$sinfo %>% select(unique_sample_name, vacc_inf), by = "unique_sample_name") %>%
    mutate(vacc_inf = factor(vacc_inf, levels = c("V- I-", "V+ I-", "V- I+", "V+ I+")))
  
  plt_out <- ggplot(plt_dat, aes(x = V1, y = V2, colour = vacc_inf)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_colour_manual(values = RColorBrewer::brewer.pal(9, "Purples")[seq(3, 9, 2)] %>%
                          setNames(c("V- I-", "V+ I-", "V- I+", "V+ I+"))) +
    labs(colour = "Self-reported\nimmune status",
         x = "UMAP1", y = "UMAP2") +
    theme_classic(16)
  
  return(plt_out)
}

# Protein IQR plot
make_prot_iqr_plot <- function(prot_in) {
  plt_dat <- prot_in$npm %>%
    pivot_longer(cols = -unique_sample_name) %>%
    summarise(iqr = IQR(value), .by = name) %>%
    # Order by IQR
    arrange(iqr) %>%
    mutate(name = factor(name, unique(name)))
  
  plt_out <- ggplot(plt_dat, aes(x = iqr)) +
    geom_histogram(bins = 30, colour = "black", fill = "grey") +
    labs(x = "IQR", y = "Protein count (N)") +
    theme_classic(16) +
    theme()
  
  return(list("data" = plt_dat, "plot" = plt_out))
}

# Combine QC plots into one
make_prot_qc_plot <- function(outliers_a, outliers_o, prot_det, prot_iqr, prot_cv, prot_pca) {
  dsgn <- "
  AD
  BE
  CF
  "
  
  plt_out <- wrap_plots(
    outliers_a$plots[[1]] + theme_classic(16),
    outliers_o$plots[[1]] + theme_classic(16),
    prot_pca,
    prot_det$plot,
    prot_cv$plot,
    prot_iqr$plot,
    design = dsgn
  ) + plot_layout(axis_titles = "collect")
  
  return(plt_out)
}

# Protein associations with clinical variables
run_prot_clin_assoc <- function(prot_in, ukb_assoc1) {
  clinvar <- c("age_group", "sex")
  
  test_results <- prot_in$npm %>%
    left_join(prot_in$sinfo %>% select(unique_sample_name, !!clinvar),
              by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, !!clinvar),
                 names_to = "protein", values_to = "npm") %>%
    # Treat age and sex as numeric variables. This makes no difference for sex as it contains only two categories. For age_group the variable is ordinal and evenly spaced
    mutate(across(!!clinvar, \(x) as.numeric(factor(x, levels = sort(unique(x)))))) %>%
    pivot_longer(cols = -c(unique_sample_name, protein, npm),
                 names_to = "clinvar", values_to = "value") %>%
    group_by(protein, clinvar) %>% nest() %>%
    # Run model on each protein
    mutate(lm_res = map(data, \(x) {
      lm(npm ~ value, data = x) %>%
        summary() %>%
        pluck("coefficients") %>%
        as.data.frame() %>%
        slice(2)
      
    })) %>%
    # Extract estimates and p-values
    mutate(estimate = map(lm_res, "Estimate"),
           p_value = map(lm_res, "Pr(>|t|)")) %>%
    unnest(c(estimate, p_value)) %>%
    ungroup() %>% select(-data, -lm_res) %>%
    mutate(fdr = p.adjust(p_value, method = "fdr"), .by = clinvar)
  
  # Add info on age associations found in UKB
  test_results <- test_results %>%
    left_join(prot_in$binfo %>% select(unique_protein_name, uniprot),
              by = c("protein" = "unique_protein_name")) %>%
    left_join(ukb_assoc1 %>% select(uniprot, age_logp, sex_logp) %>%
                pivot_longer(cols = -uniprot, names_to = "clinvar", values_to = "ukb_logp") %>%
                mutate(clinvar = case_when(clinvar == "age_logp" ~ "age_group",
                                           clinvar == "sex_logp" ~ "sex")),
              by = c("uniprot", "clinvar"), relationship = "many-to-many")

  return(test_results)
}

# Save a table of the associations
save_prot_clin_assoc <- function(prot_clin_assoc, prot_seroclust_assoc, ifn_prot_hms,
                                 filename = "prot_clin_assoc.csv") {
  
  # Mark proteins that pop up in protein-serocluster or protein-Ab associations
  prot_clin_assoc <- mutate(
    prot_clin_assoc,
    assoc_serocluster = case_when(protein %in% (prot_seroclust_assoc$kruskal_wallis %>%
                                      filter(fdr < 0.05) %>% pull(protein)) ~ 1),
    concordant_assoc_abs = case_when(protein %in% (ifn_prot_hms$proteins_p %>%
                                       pull(protein) %>% unique() %>%
                                       rep(., each = 2) %>% paste0(c("_1", "_2"))) ~ 1)
  )
  
  file_name <- make_fn("protein_clinical_association", filename)
  write_csv(prot_clin_assoc, file_name)
  
  # Return the file name
  return(file_name)
}

# Plot protein-clinical associations
plot_prot_clin_assoc <- function(prot_clin_assoc_tbl) {
  plt_dat <- read_csv(prot_clin_assoc_tbl, show_col_types = F) %>%
    mutate(clinvar = clinvar %>% str_replace_all("_", " ") %>% str_to_sentence())
  
  plt <- plt_dat %>%
    # Combine annotation to make one variable for colouring
    mutate(assoc = case_when(assoc_serocluster == 1 & concordant_assoc_abs == 1 ~ "Both", assoc_serocluster == 1 ~ "Serocluster", concordant_assoc_abs == 1 ~ "Abs", T ~ "None")) %>%
    mutate(assoc = factor(assoc, levels = c("None", "Abs", "Serocluster", "Both"))) %>%
    # Sort to change drawing order, to avoid obscuring positive points with negative ones
    arrange(assoc) %>%
    ggplot(aes(x = estimate, y = -log10(fdr))) +
    geom_point(aes(colour = assoc), size = 2, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
    geom_text_repel(aes(label = protein), data = plt_dat %>% slice_min(fdr, n = 10, by = clinvar) %>%
                      mutate(protein = str_remove(protein, "_\\d$")), max.overlaps = 25) +
    scale_colour_manual(values = c("None" = "lavenderblush2", "Both" = "gold2", "Serocluster" = "sienna2", "Abs" = "lightskyblue")) +
    scale_x_continuous(limits = center_limits()) +
    facet_wrap(~ clinvar) +
    labs(colour = "Associations", x = "Beta") +
    theme_classic(14)
  
  return(plt)
}

# Make plots comparing betas from protein-questionnaire associations of the two regions
make_clin_beta_plot <- function(clin_assoc_region) {
  # Bind together the results from the two regions for plotting
  plt_dat <- rbind(clin_assoc_region$Stockholm %>% mutate(region = "Stockholm"),
                   clin_assoc_region$Gothenburg %>% mutate(region = "Gothenburg")) %>%
    select(protein, clinvar, estimate, p_value, region) %>%
    distinct()
  
  plt1 <- plt_dat %>%
    # Make wider data to put one region on each axis
    pivot_wider(id_cols = c(protein, clinvar), names_from = region, values_from = estimate) %>%
    ggplot(aes(x = Stockholm, y = Gothenburg)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    scale_x_continuous(limits = ggh4x::center_limits()) +
    scale_y_continuous(limits = ggh4x::center_limits()) +
    facet_wrap(~ clinvar, scales = "free") +
    labs(x = "Beta Stockholm", y = "Beta Gothenburg") +
    theme_classic(14)
  
  plt2 <- plt_dat %>%
    # Include only proteins with a p-value below 0.05 in both regions
    filter(p_value < 0.05) %>%
    pivot_wider(id_cols = c(protein, clinvar), names_from = region, values_from = estimate) %>%
    filter(!is.na(Stockholm), !is.na(Gothenburg)) %>%
    ggplot(aes(x = Stockholm, y = Gothenburg)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    scale_x_continuous(limits = ggh4x::center_limits()) +
    scale_y_continuous(limits = ggh4x::center_limits()) +
    facet_wrap(~ clinvar, scales = "free") +
    labs(x = "Beta Stockholm", y = "Beta Gothenburg") +
    theme_classic(14)
  
  return(list("all" = plt1, "signif" = plt2))
}

# Proteins and serology measurements associations (linear regression)
run_prot_ab_assoc <- function(prot_in, serol_seropos) {
  
  test_results <- serol_seropos$per_ab %>%
    mutate(name = rename_serol(name)) %>%
    select(unique_sample_name, name, seroclass) %>%
    left_join(prot_in$npm, by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, name, seroclass),
                 names_to = "protein", values_to = "npm") %>%
    filter(!is.na(npm)) %>%
    mutate(seroclass = as.factor(seroclass)) %>%
    group_by(name, protein) %>% nest() %>%
    mutate(p_value = map(data, \(x) {
      glm(seroclass ~ npm, family = "binomial", data = x) %>%
        summary() %>% pluck("coefficients") %>% as_tibble() %>%
        pluck("Pr(>|z|)", 2)
    })) %>%
    # Add an estimate of direction equal to the difference in median NPM between the groups, positive if the seropositive group is higher than the seronegative
    mutate(estimate = map(data, \(x) {
      median(x %>% filter(seroclass == "Positive") %>% pull(npm)) -
        median(x %>% filter(seroclass == "Negative") %>% pull(npm))
    })) %>%
    select(-data) %>%
    unnest(c(estimate, p_value)) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr"), .by = name) %>%
    select(protein, ab = name, estimate, p_value, fdr) %>%
    arrange(protein, ab)

  # Make heatmap of results
  # Split it into multiple parts for visibility
  splits <- ceiling(seq_along(unique(test_results$protein)) /
                      ceiling(length(unique(test_results$protein)) / 5)) %>%
    setNames(unique(test_results$protein))

  plt_dat <- test_results %>%
    mutate(fill_p = case_when(p_value < 0.05 & estimate > 0 ~ "Positive",
                              p_value < 0.05 & estimate < 0 ~ "Negative",
                              p_value > 0.05 ~ "None"),
           fill_fdr = case_when(fdr < 0.05 & estimate > 0 ~ "Positive",
                                fdr < 0.05 & estimate < 0 ~ "Negative",
                                fdr > 0.05 ~ "None")) %>%
    select(protein, ab, fill_p, fill_fdr, p_value, fdr)

  test_hm_p <- plt_dat %>%
    mutate(facet = splits[protein]) %>%
    # Reorder proteins to start at the top
    arrange(facet, desc(protein)) %>%
    mutate(protein = factor(protein, levels = unique(protein))) %>%
    select(protein, ab, fill_p, facet) %>%
    ggplot(aes(x = ab, y = protein, fill = fill_p)) +
    geom_tile(colour = "grey", linewidth = 0.1) +
    facet_grid2(cols = vars(facet), scales = "free", drop = T, independent = "all") +
    labs(x = NULL, y = NULL, fill = "Association") +
    scale_fill_manual(values = c("Positive" = "sienna2", "Negative" = "skyblue2", "None" = "white")) +
    theme_minimal() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  test_hm_fdr <- plt_dat %>%
    mutate(facet = splits[protein]) %>%
    # Reorder proteins to start at the top
    arrange(facet, desc(protein)) %>%
    mutate(protein = factor(protein, levels = unique(protein))) %>%
    select(protein, ab, fill_fdr, facet) %>%
    ggplot(aes(x = ab, y = protein, fill = fill_fdr)) +
    geom_tile(colour = "grey", linewidth = 0.1) +
    facet_grid2(cols = vars(facet), scales = "free", drop = T, independent = "all") +
    labs(x = NULL, y = NULL, fill = "Association") +
    scale_fill_manual(values = c("Positive" = "sienna2", "Negative" = "skyblue2", "None" = "white")) +
    theme_minimal() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list("assoc" = test_results,
              "heatmap_p" = test_hm_p,
              "heatmap_fdr" = test_hm_fdr))
}

# Proteins and serology cluster associations (kruskal-wallis)
run_prot_seroclust_assoc <- function(prot_in, seroclust) {
  # Kruskal-Wallis
  kwp <- prot_in$npm %>%
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster),
              by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, serocluster),
                 names_to = "protein", values_to = "npm") %>%
    group_by(protein) %>% nest() %>%
    mutate(p_value = map(data, \(x) {
      kruskal.test(npm ~ serocluster, data = x)$p.value
    })) %>%
    # select(-data) %>%
    unnest(p_value) %>%
    # ungroup(protein) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    left_join(prot_in$binfo %>% select(unique_protein_name, assay_method = method),
              by = c("protein" = "unique_protein_name"))
  
  # Wilcoxon
  wp <- kwp %>%
    filter(fdr < 0.05)
  
  if (nrow(wp) > 0) {
    wp <- wp %>%
      select(-p_value, -fdr) %>%
      group_by(protein) %>%
      mutate(p_value = map(data, \(x) {
        test_results <- pairwise.wilcox.test(x$npm, x$serocluster, p.adjust.method = "none")$p.value %>%
          as.data.frame() %>%
          # grp1 from columns, grp2 from rows, so that the order makes more sense (lower numbers in grp 1)
          rownames_to_column("grp2") %>%
          pivot_longer(cols = -grp2, names_to = "grp1", values_to = "p_value") %>%
          filter(!is.na(p_value)) %>%
          relocate(grp1, grp2) %>%
          mutate(protein = protein, .before = 1)
        
        # Compute differences in median signal between seroclusters, add to output
        grp_med <- summarise(x, med = median(npm), .by = serocluster)
        grp_diff <- grp_med %>% select(grp1 = serocluster, grp2 = serocluster) %>%
          complete(grp1, grp2) %>%
          filter(grp1 != grp2) %>%
          mutate(med_diff = grp_med[match(grp1, grp_med$serocluster), "med", drop = T] -
                   grp_med[match(grp2, grp_med$serocluster), "med", drop = T])
        
        return(test_results %>% left_join(grp_diff, by = c("grp1", "grp2")))
      })) %>%
      pull(p_value) %>%
      bind_rows() %>%
      mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
      left_join(prot_in$binfo %>% select(unique_protein_name, assay_method = method),
                by = c("protein" = "unique_protein_name")) %>%
      ungroup()
    
  }
  
  kw_qq <- kwp %>%
    select(protein, p_value) %>%
    # Observed and expected p-values for plot
    mutate(obs = -log10(sort(p_value)),
           exp = -log10(1:length(p_value) / (length(p_value) + 1))) %>%
    ggplot(aes(x = exp, y = obs)) +
    geom_point(shape = 1) +
    geom_abline(slope = 1) +
    scale_colour_brewer(palette = "Set2") +
    theme_classic() +
    theme(aspect.ratio = 1)
  
  return(list("kruskal_wallis" = kwp %>% ungroup() %>% select(-data),
              "wilcoxon" = wp))
}

# Make plots for the protein serocluster associations
make_prot_seroclust_plots <- function(prot_seroclust_assoc, prot_dat, seroclust, clust_cols) {
  # QQ plots of Kruskal-Wallis p-values
  kw_qq <- prot_seroclust_assoc$kruskal_wallis %>%
    select(protein, p_value) %>%
    # Observed and expected p-values for plot
    mutate(obs = -log10(sort(p_value)),
           exp = -log10(1:length(p_value) / (length(p_value) + 1))) %>%
    ggplot(aes(x = exp, y = obs)) +
    geom_point(shape = 1) +
    geom_abline(slope = 1) +
    scale_colour_brewer(palette = "Set2") +
    theme_classic() +
    theme(aspect.ratio = 1)
  
  # Heatmap and boxplots of Wilcoxon p-values
  if (nrow(prot_seroclust_assoc$wilcoxon) > 0) {
    # Heatmap for wilcoxon test results comparing methods
    wilcox_hm <- prot_seroclust_assoc$wilcoxon %>%
      mutate(significant = fdr < 0.05) %>%
      # If FDR significant colour by median difference, otherwise NA
      mutate(med_diff = case_when(significant ~ med_diff)) %>%
      unite("grp_pair", grp1, grp2, sep = "_") %>%
      select(protein, grp_pair, med_diff, significant) %>%
      filter(!is.na(med_diff)) %>%
      ggplot(aes(x = grp_pair, y = protein, fill = med_diff)) +
      geom_tile() +
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "RdBu"))) +
      theme_minimal()
    
    # Boxplots of proteins with kruskal-wallis fdr < 0.05
    box_prots <- prot_seroclust_assoc$kruskal_wallis %>% filter(fdr < 0.05) %>% pull(protein)
    wilcox_box <- prot_dat$npm %>%
      select(unique_sample_name, !!box_prots) %>%
      left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster), by = "unique_sample_name") %>%
      pivot_longer(cols = -c(unique_sample_name, serocluster),
                   names_to = "protein", values_to = "npm") %>%
      group_by(protein) %>% nest() %>%
      mutate(boxplt = map(data, \(x) {
        # Current protein, to avoid name clashes
        prot <- protein
        
        # FDR values for labels, convert to stars for a less crowded plot
        fdr_df <- prot_seroclust_assoc$wilcoxon %>% filter(protein == prot, fdr < 0.05)
        fdr_comparisons <- apply(fdr_df, 1, \(rw) {c(rw["grp1"], rw["grp2"])}, simplify = F)
        
        fdr_labs <- paste0("p == ", signif(fdr_df$fdr, 2))
        
        ggplot(x, aes(x = serocluster, y = npm)) %>%
          geom_beebox(aes_bee = aes(colour = serocluster), param_bee = list(cex = 1.5),
                      param_box = list(colour = "black", width = 0.25, linewidth = 0.7,
                                       alpha = 0, show.legend = F)) +
          geom_signif(comparisons = fdr_comparisons, annotations = fdr_labs,
                      step_increase = 0.25, colour = "black", parse = T) +
          scale_colour_manual(values = clust_cols %>% setNames(1:4)) +
          labs(x = NULL, y = "Protein levels [AU]", title = str_remove(prot, "_\\d$"), colour = "Serocluster") +
          # Prevent clipping of the p-values when combining the plots
          coord_cartesian(clip = "off") +
          # Increase point size in legend
          guides(colour = guide_legend(override.aes = list(size = 5))) +
          theme_classic(13) +
          theme(plot.title = element_text(hjust = 0, size = 12), axis.text = element_text(size = 12),
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                axis.title = element_text(size = 12),
                legend.direction = "horizontal", legend.position = "bottom",
                legend.title = element_text(size = 12), legend.text = element_text(size = 12))
      })) %>%
      select(protein, boxplt) %>%
      arrange(protein)
    
  } else {
    wilcox_hm <- "No Kruskal-Wallis FDR values < 0.05"
    wilcox_box <- "No Kruskal-Wallis FDR values < 0.05"
  }
  
  
  return(list("qq_plots" = kw_qq,
              "wilcox_hm" = wilcox_hm,
              "wilcox_box" = wilcox_box))
}

# Save the boxplots as a pdf
save_prot_seroclust_boxplts <- function(prot_seroclust_plts) {
  filename <- make_fn("serocluster_comparison", "seroclust_vs_protein.pdf")
  
  # Combine to make plot
  plts <- prot_seroclust_plts$wilcox_box$boxplt
  plt <- lapply(split(seq_along(plts), rep(1:(length(plts) / 3), each = 3)), \(i) {
    wrap_plots(plts[i], nrow = 1) + plot_layout(axis_titles = "collect") & theme(legend.position = "bottom")
  }) %>% wrap_plots(ncol = 1, guides = "collect") & theme(legend.position = "bottom")
  
  ggsave(filename, width = 7, height = 12, plot = plt)
  
  # File name as output to set format = "file" in the target
  return(filename)
}

# Machine learning (Regularised regression (Lasso) or gradient boosted trees (XGBoost)) for predicting seroclusters using circulating proteins
run_prot_ml <- function(prot_dat, seroclust, train_test_split, clust_cols, engine = "glmnet", method = "lasso", tune_v = 5, tune_grid_size = 30, rseed = 123) {
  # train_test_split, numeric of length 1 specifying the proportion of training samples. Alternatively a character vector (e.g. "stockholm", "gothenburg") in which case the split is done simply by taking one region as training and the other as testing set (matches with region to use as training set).
  # clust_cols, colour to use for clusters in ROC curves, same length as number of classes
  
  npm_in <- prot_dat$npm %>%
    # Add cluster memberships (the outcome)
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster),
              by = "unique_sample_name") %>%
    mutate(serocluster = as.factor(serocluster))
  
  # Get indices of training samples if splitting by region
  if (is.character(train_test_split)) {
    split_in <- match(prot_dat$sinfo %>%
                        filter(str_detect(region, regex(paste0("^", train_test_split), ignore_case = T))) %>%
                        pull(unique_sample_name),
                      npm_in$unique_sample_name)
    
  } else {
    split_in <- train_test_split
  }
  
  ml_results <- tidy_classification(x = npm_in,
                                    response = serocluster,
                                    train_test_split = split_in,
                                    engine = engine, method = method,
                                    preproc_steps = . %>%
                                      update_role(unique_sample_name, new_role = "ID") %>%
                                      step_zv(all_predictors()) %>%
                                      step_corr(all_predictors(), threshold = 0.6) %>%
                                      step_normalize(all_predictors()),
                                    tune_v = tune_v, tune_grid_size = tune_grid_size,
                                    rseed = 123)
  # Change colours to match other plots
  ml_results$roc_curve <- ml_results$roc_curve +
    scale_colour_manual(values = clust_cols)
  
  return(ml_results) 
}

# Also ML for predicting sex, age, and region (questionnaire variables)
run_prot_ml_others <- function(prot_dat, train_test_split, engine = "glmnet", method = "lasso", tune_v = 5, tune_grid_size = 30, rseed = 123) {
  q_vars <- c("sex", "age_group", "region", "vaccinated", "infected")
  dat_in <- prot_dat$npm %>%
    left_join(prot_dat$sinfo %>% select(unique_sample_name, !!q_vars), by = "unique_sample_name") %>%
    mutate(vaccinated = case_when(vaccinated > 0 ~ 1, vaccinated == 0 ~ 0)) %>%
    column_to_rownames("unique_sample_name")
  
  lassos <- sapply(q_vars, \(qvar) {
    # Order response
    dat_in[, qvar] <- factor(dat_in[[qvar]], levels = sort(unique(as.character(dat_in[[qvar]]))))
    # Omit NA
    dat_in %>% select(everything(), -!!q_vars, !!qvar) %>%
      filter(complete.cases(.)) %>%
      tidy_classification(response = !!qvar,
                          train_test_split = train_test_split, engine = engine, method = method,
                          tune_v = tune_v, tune_grid_size = tune_grid_size, rseed = rseed)
  }, simplify = F, USE.NAMES = T)
  
  return(lassos)
}

# Also ML for predicting vaccination on age-adjusted protein levels
run_prot_ml_vacc <- function(prot_dat, train_test_split = 0.7, engine = "glmnet", method = "lasso", tune_v = 5, tune_grid_size = 30, rseed = 123) {
  dat_in <- prot_dat$npm %>%
    # Add sample info
    left_join(prot_dat$sinfo %>% select(unique_sample_name, age_group, vaccinated), by = "unique_sample_name") %>%
    filter(!is.na(vaccinated)) %>%
    pivot_longer(cols = -c(unique_sample_name, age_group, vaccinated)) %>%
    # Adjust each protein for age using linear model 
    nest_by(name) %>%
    mutate(value_adj = list(
      lm(value ~ age_group, data = data) %>% resid()
    )) %>%
    unnest(c(data, value_adj)) %>%
    # Convert to wide format again and change vaccination to binary factor
    pivot_wider(id_cols = c("unique_sample_name", "vaccinated"), names_from = "name", values_from = "value_adj") %>%
    mutate(vaccinated = case_when(vaccinated > 0 ~ "1", vaccinated == 0 ~ "0"),
           vaccinated = as.factor(vaccinated)) %>%
    column_to_rownames("unique_sample_name")
  
  lasso <- tidy_classification(x = dat_in, response = vaccinated,
                               train_test_split = train_test_split,
                               engine = engine, method = method, tune_v = tune_v,
                               tune_grid_size, rseed = rseed)
  return(lasso)
}

#  CRP classes
classify_crp <- function(prot_dat, n_sd = 3) {
  # n_sd, number of SDs to use for cutoff (provide multiple values for multiple classes)
  
  dat <- prot_dat$npm %>%
    select(unique_sample_name, CRP)
  # Estaimted negative proportion
  np <- prot_dat$sinfo %>%
    mutate(np = 1 - (sum(infected) / n())) %>%
    pull(np) %>%
    unique()
  
  # Get classifications and cutoff
  cls <- lapply(n_sd, \(nsd) {
    dens_results <- density_cutoff(x = dat$CRP, n_sd = nsd, neg_prop = np)
    
    return(data.frame(unique_sample_name = dat$unique_sample_name,
                      CRP = dat$CRP,
                      class = dens_results$class,
                      n_sd = nsd,
                      cutoff = dens_results$cutoff))
  }) %>% bind_rows() %>%
    # Get highest cutoff wherer a sample is positive
    summarise(CRP = unique(CRP),
              cutoff = max(c(-Inf, cutoff[class])),
              class = case_when(!any(class) ~ "Negative",
                                # Include -Inf to avoid warning where max returns -Inf when there are no TRUE values (no value above any of the cutoffs)
                                T ~ paste0(">", max(c(-Inf, n_sd[class])), " SD")),
              .by = unique_sample_name)
  
  return(cls)
}

# Plot the CRP classification
make_crp_plot <- function(crp_classes, serol_bbnorm, palette = "Oranges") {
  # Beeswarm plot of classification
  crp_plt <- crp_classes %>%
    ggplot(aes(x = 0, y = CRP)) +
    geom_beeswarm(aes(colour = class), size = 1.5) +
    geom_hline(aes(yintercept = cutoff), colour = "grey", linetype = "dashed") +
    geom_boxplot(alpha = 0, outliers = F, width = 0.025, linewidth = 1) +
    scale_colour_brewer(palette = palette) +
    labs(y = "CRP level [AU]", colour = "Class") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  igm_plt <- serol_bbnorm$mfi %>%
    select(unique_sample_name, IgM = contains("human_igm")) %>%
    left_join(crp_classes %>% select(unique_sample_name, class), by = "unique_sample_name") %>%
    filter(!is.na(class)) %>%
    ggplot(aes(x = 0, y = IgM)) +
    geom_beeswarm(aes(colour = class), size = 1.5) +
    geom_boxplot(alpha = 0, outliers = F, width = 0.025, linewidth = 1) +
    labs(colour = "CRP class") +
    scale_colour_brewer(palette = palette) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  filename <- make_fn("crp_analysis", "crp_plots.pdf")
  pdf(file = filename, onefile = T, height = 2.5, width = 4)
  print(crp_plt)
  print(igm_plt)
  dev.off()
  
  return(filename)
}

characterise_crp <- function(crp_classes, prot_dat, seroclust, prot_clust) {
  plt_dat <- crp_classes %>%
    left_join(prot_dat$sinfo %>%
                select(unique_sample_name, age_group, sex, region, infected, vaccinated, inf_diff, vacc_diff, vacc_inf, vaccine_company),
              by = "unique_sample_name") %>%
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster), by = "unique_sample_name") %>%
    mutate(infected = as.character(infected), vaccinated = as.character(vaccinated)) %>%
    left_join(prot_clust$clust_dat, by = "unique_sample_name") %>%
    rename(proteocluster = cluster) %>%
    mutate(proteocluster = as.character(proteocluster))
  
  plt <- plt_dat %>% select(-c(unique_sample_name, CRP, cutoff, class)) %>%
    colnames() %>%
    sapply(\(x) {
      comp_groups(as.formula(paste0(x, " ~ class")), plt_dat, sim_p = T, discr_palette = case_when(x == "proteocluster" ~ "Oranges", T ~ "Purples"))$plt +
        labs(x = "CRP class",
             # Better names for numerical variables
             y = case_when(x == "inf_diff" ~ "Months from infection",
                           x == "vacc_diff" ~ "Months from vaccination",
                           T ~ "Proportion"),
             # Better names for categorical variables
             fill = case_when(x == "vacc_inf" ~ "Self-reported immune status",
                              T ~ x %>% str_replace("_", " ") %>% str_to_sentence()))
    }, simplify = F)
  
  filename <- make_fn("crp_analysis", "crp_characteristics.pdf")
  ggsave(filename = filename, plot = wrap_plots(plt), height = 10, width = 16)
  
  return(filename)
}

# Run Weighted Gene Correlation Network Analysis (WGCNA) to get protein modules
run_wgcna <- function(prot_dat) {
  prot_in <- prot_dat$npm %>% column_to_rownames("unique_sample_name")
  
  power_check <- pickSoftThreshold(data = prot_in, dataIsExpr = T, verbose = 0)
  
  power_plts <- power_check$fitIndices %>% rename_with(tolower) %>%
    dplyr::select(1, 2, 5:7) %>%
    pivot_longer(cols = -power) %>%
    ggplot(aes(x = power, y = value)) +
    geom_point() + geom_line(aes(group = name)) +
    geom_vline(xintercept = power_check$powerEstimate, colour = "red", linetype = "dotted") +
    facet_wrap(~ name, scales = "free_y") +
    theme_classic()
  
  wgcna <- blockwiseModules(datExpr = prot_in, corType = "bicor",
                            power = power_check$powerEstimate,
                            verbose = 0, saveTOMs = F)
  
  # Heatmap of the module 
  me_heatmap <- Heatmap(wgcna$MEs)
  
  return(list("power" = power_check$powerEstimate, "power_plots" = power_plts, "wgcna" = wgcna))
}

# Visualise the protein WGCNA modules
visualise_wgcna_modules <- function(prot_dat, prot_wgcna, seroclust) {
  
  # Heatmap of proteins and their adjacency
  # Make topology overlap matrix
  tom <- TOMsimilarityFromExpr(prot_dat$npm %>% column_to_rownames("unique_sample_name"),
                               power = prot_wgcna$power)
  # Get dissimilarity and transform for plotting
  plot_tom <- (1 - tom) ** 10
  
  diag(plot_tom) <- NA
  
  # Draw with pheatmap to assign to a variable
  # Use modules and their colours to annotate columns
  annot_df <- data.frame(module = prot_wgcna$wgcna$colors, row.names = names(prot_wgcna$wgcna$colors))
  annot_col <- list(module = unique(prot_wgcna$wgcna$colors) %>% setNames(unique(prot_wgcna$wgcna$colors)))
  tom_hm <- pheatmap(plot_tom, color = rev(hcl.colors(100)), show_rownames = F, show_colnames = F,
                     cluster_rows = prot_wgcna$wgcna$dendrograms[[1]],
                     cluster_cols = prot_wgcna$wgcna$dendrograms[[1]],
                     treeheight_row = 0, annotation_col = annot_df, annotation_colors = annot_col,
                     name = "Dissimilarity")
  
  # Module eigengene (MEs) adjacency (defined as (1 + cor) / 2)
  me_adj <- prot_wgcna$wgcna$MEs %>%
    cor() %>% add(1) %>% divide_by(2) %>%
    pheatmap(color = hcl.colors(10), breaks = seq(0, 1, length.out = 10),
             border_color = "white", display_numbers = T, fontsize_number = 12,
             name = "Adjacency")
  
  return(list("tom_hm" = tom_hm, "me_adj" = me_adj))
}

# Perform PCA per module to get things like variable contribution and PC importance, plot differences in proteotypes per module
wgcna_module_pca <- function(prot_dat, prot_wgcna, prot_clust_effect) {
  sapply(unique(prot_wgcna$wgcna$colors), \(x) {
    prots <- names(prot_wgcna$wgcna$colors[prot_wgcna$wgcna$colors == x])
    pca <- prcomp(prot_dat$npm %>% select(!!prots),
                  center = T, scale. = T) %>%
      summary()
    
    # Make a table of all proteins and their (rank of) importance for PC1
    pca$prot_imp <- pca$rotation %>% as.data.frame() %>%
      select(PC1) %>% rownames_to_column("unique_protein_name") %>%
      # Add protein name and uniprot id
      mutate(protein = str_remove(unique_protein_name, "_\\d$")) %>%
      left_join(prot_dat$binfo %>% select(unique_protein_name, uniprot),
                by = "unique_protein_name") %>%
      # Add rank
      mutate(importance_rank = rank(-abs(PC1))) %>%
      arrange(importance_rank) %>%
      as_tibble()
    
    # Pick out top 10 most contributing proteins to PC1
    pca$top_prot <- pca$prot_imp %>%
      slice_min(order_by = importance_rank, n = 10)
    
    pca$effect_plot <- make_prot_clust_effect_plots(prot_clust_effect %>%
                                                      filter(name %in% pca$top_prot$unique_protein_name) %>%
                                                      # Order in order of importance
                                                      mutate(name = factor(name, levels = rev(pca$top_prot$unique_protein_name))) %>%
                                                      mutate(comparison = factor(comparison, levels = as.character(sort(unique(as.numeric(comparison)))))),
                                                    spread = 0.6, point_size = 0.5)
    
    return(pca)
  }, simplify = F)
}

# Check how proteins overlapping between technologies end up in WGCNA
# Code provided by Maria Bueno Alvez (modified to add figure labels and move the legend)
check_wgcna_overlap <- function(module_proteins) {
  # Theme & palette
  theme_hpa <- 
    function(angled = F, axis_x = T, axis_y = T, facet_title = T) {
      t <- 
        theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.2, "lines"), 
          panel.background=element_rect(fill="white"),
          panel.border = element_blank(),
          plot.title = element_text(face = "bold",
                                    size = rel(1), hjust = 0.5),
          plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.ticks.length = unit(.25, "cm"),
          axis.line = element_line(linewidth = 0.5),
          axis.text = element_text(size = rel(1), color = 'black'),
          legend.key = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=rel(0.8)),
          legend.key.size= unit(0.7, "cm"),
          legend.title = element_text(size=rel(1)),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="grey90",fill="grey90"),
          strip.text = element_text(face="bold")
        )
      
      if(angled) {
        t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
      }
      
      if(axis_x == F) {
        t <- t +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                axis.title.x = element_blank())
      } 
      
      if(axis_y == F) {
        t <- t +
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                axis.title.y = element_blank())
      }
      if(facet_title == F) {
        t <- t + theme(strip.text = element_blank())
      }
      return(t)
    }
  
  # Upper case first letter in modules
  module_proteins$module <- str_to_sentence(module_proteins$module)
  
  pal_module <- c(
    Grey = "#D9D9D9",
    Yellow = "#EEED36",
    Turquoise = "#7ADDD0",
    Green = "#ACCB62",
    Brown = "#C18747",
    Blue = "#5377FA"
  )
  
  
  # Supplementary plot
  p1 <- 
    module_proteins |> 
    count(module) |> 
    mutate(module = factor(module, levels = names(pal_module))) |> 
    ggplot(aes(module, n, fill = module)) +
    geom_col(show.legend = F) +
    geom_text(aes(label = n), vjust = -0.5) +
    scale_fill_manual(values = pal_module) +
    theme_hpa(axis_x = F) +
    ylab("Number of assays") +
    coord_cartesian(clip = "off") +
    theme(axis.x.ticks = element_blank(),
          axis.x.text = element_blank(),
          axis.line.x = element_line(),
          axis.title.x = element_text()) +
    xlab("Module") 
  
  p2 <- 
    module_proteins |> 
    distinct(assay, module) |> 
    count(module) |> 
    mutate(module = factor(module, levels = names(pal_module))) |> 
    ggplot(aes(module, n, fill = module)) +
    geom_col(show.legend = F) +
    geom_text(aes(label = n), vjust = -0.5) +
    scale_fill_manual(values = pal_module) +
    theme_hpa(axis_x = F) +
    ylab("Number of unique assays") +
    coord_cartesian(clip = "off") +
    theme(axis.x.ticks = element_blank(),
          axis.x.text = element_blank(),
          axis.line.x = element_line(),
          axis.title.x = element_text()) +
    xlab("Module") 
  
  # Identify common assays Olink & Alamar
  duplicated_assays <- 
    module_proteins |> 
    group_by(assay) |> 
    count(assay) |> 
    filter(n > 1) 
  
  # Step 1: Flag assays with consistent modules across both technologies
  module_ordered <- 
    module_proteins |>
    filter(assay %in% duplicated_assays$assay) |>
    mutate(Technology = sub(".*_", "", unique_assay_name),
           module = factor(module, levels = rev(names(pal_module))))|>
    group_by(assay, module) |>
    mutate(same_module_across_technologies = n_distinct(Technology) > 1) |>
    ungroup() |>
    arrange(desc(same_module_across_technologies), module, assay)
  
  # Step 2: Plot with updated ordering
  p3 <- 
    module_ordered |>
    ggplot(aes(y = factor(assay, levels = rev(unique(assay))), 
               x = Technology, fill = module)) +
    geom_tile() +
    scale_fill_manual(values = pal_module) +
    coord_fixed() +
    theme_hpa(angled = TRUE, axis_x = F) +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    guides(fill = guide_legend(ncol = 2)) +
    labs(x = NULL, y = NULL, fill = "Module")
  
  # Combine
  dsgn <- "aad\nbbd\nccd"
  plt <- wrap_plots(p1 + labs(tag = "a"),
                    p2 + labs(tag = "b"),
                    guide_area(),
                    p3 + labs(tag = "c"),
                    design = dsgn,
                    guides = "collect")
  
  return(list("data" = module_ordered, "plt" = plt))
}

# Check associations of WGCNA module proteins with UKB traits
make_wgcna_ukb_table <- function(prot_dat, prot_wgcna, ukb1, ukb2) {
  wgcna_ukb <- sapply(unique(prot_wgcna$wgcna$colors), \(x) {
    prots <- names(prot_wgcna$wgcna$colors[prot_wgcna$wgcna$colors == x])
    tot_df <- ukb1 %>% select(assay, olinkid, matches("beta|logp")) %>%
      full_join(ukb2 %>% select(assay, olinkid, matches("beta|logp")),
                relationship = "many-to-many", by = c("assay", "olinkid")) %>%
      # Left join to look only at proteins available in the UKB study
      left_join(prot_dat$binfo %>% select(unique_protein_name, assay),
                relationship = "many-to-many", by = "assay") %>%
      filter(unique_protein_name %in% prots) %>%
      relocate(unique_protein_name)
    # To simplify the code, divide into two, convert to long format separately and merge again
    beta_df <- tot_df %>%
      select(unique_protein_name, assay, olinkid, contains("_beta")) %>%
      pivot_longer(cols = -c(unique_protein_name, assay, olinkid), names_to = "trait", values_to = "beta") %>%
      mutate(trait = str_remove(trait, "_beta"))
    logp_df <- tot_df %>%
      select(unique_protein_name, assay, olinkid, contains("_logp")) %>%
      pivot_longer(cols = -c(unique_protein_name, assay, olinkid), names_to = "trait", values_to = "logp") %>%
      mutate(trait = str_remove(trait, "_logp"))
    out_df <- beta_df %>%
      left_join(logp_df, by = c("unique_protein_name", "assay", "olinkid", "trait"),
                relationship = "many-to-many") %>%
      # Number the variables for nicer axes in plot
      mutate(trait_num = as.integer(factor(trait, levels = unique(trait)))) %>%
      distinct()
    
    return(out_df)
  }, simplify = F)
  
  return(wgcna_ukb)
}

# Make plots visualising associations from UKB
make_wgcna_ukb_plots <- function(wgcna_ukb) {
  bar_plt <- imap(wgcna_ukb, \(x, nm) {
    x %>%
      # Some proteins were measured multiple times in the UKB study and in this study
      # Get just one p-value per protein, take the highest logp as not all proteins have overlapping OlinkIDs
      filter(logp == max(logp), .by = c(assay, trait)) %>%
      distinct(assay, trait, trait_num, logp) %>%
      count(trait_num, signif = logp > 4.769551) %>%
      mutate(prop = n / sum(n), .by = trait_num) %>%
      ggplot() +
      geom_col(aes(x = trait_num, y = prop, fill = signif), position = "stack") +
      labs(x = "Trait", y = "Proportion", title = nm, fill = "Significant") +
      scale_fill_manual(values = c("TRUE" = "#4876FF", "FALSE" = "grey")) +
      scale_x_continuous(breaks = c(1, seq(5, length(unique(x$trait_num)), by = 5))) +
      theme_classic() +
      theme(panel.grid.major.x = element_line())
  }) %>% wrap_plots(guides = "collect", axis_titles = "collect")
  
  return(bar_plt)
}

# Make plots visualising immune cell enrichment in the modules
make_wgcna_hpa_plot <- function(wgcna_module_proteins, prot_hpa_annot) {
  plt_dat <- wgcna_module_proteins %>%
    left_join(prot_hpa_annot, by = "assay") %>%
    separate_longer_delim(blood_cell_enrich, delim = ";") %>%
    count(module, blood_cell_enrich) %>%
    filter(!is.na(blood_cell_enrich)) %>%
    # Put grey module at the end of the scale
    mutate(module = factor(module, levels = c(
      append(setdiff(sort(unique(as.character(module)), decreasing = T), "grey"), "grey", after = 0)
    ))) %>%
    mutate(percent = n / sum(n) * 100, .by = blood_cell_enrich)
  
  plt_out <- ggplot(plt_dat, aes(x = percent, y = blood_cell_enrich, fill = module)) +
    geom_col() +
    scale_fill_manual(values = c("grey85", "#4876FF", "#CD853F", "#A2CD5A", "#40E0D0", "#EEEE00") %>%
                        setNames(c("grey", "blue", "brown", "green", "turquoise", "yellow"))) +
    labs(x = "Percent", y = "Blood cell type", fill = "WGCNA\nmodule") +
    theme_classic(16)
  
  return(plt_out)
}

make_wgcna_ukb_heatmap <- function(wgcna_ukb_assoc) {
  # Perform Fisher's exact test to compare number of significant associations with grey module
  tests <- imap(wgcna_ukb_assoc[!names(wgcna_ukb_assoc) == "grey"], \(x, i) {
    # Skip traits 10 and 20 that have no significant associations in any module
    x %>% filter(!trait_num %in% c(10, 20)) %>%
      # Count each protein only once within each group (main module, other modules), but if it is in both main and other count it in both
      filter(logp == max(logp), .by = c(assay, trait)) %>%
      distinct(assay, trait, trait_num, logp) %>%
      mutate(module = i) %>%
      group_by(trait, trait_num) %>% nest() %>%
      mutate(test_data = map(data, \(d) {
        tr <- trait
        # Get other module assocations
        wgcna_ukb_assoc[-which(names(wgcna_ukb_assoc) == i)] %>%
          bind_rows() %>%
          filter(trait == tr) %>%
          filter(logp == max(logp), .by = assay) %>%
          distinct(assay, trait, logp) %>%
          mutate(module = "other") %>%
          bind_rows(d)
      })) %>%
      # Count N of groups and their proportions of positive results (any associations)
      # Different traits may have different N as age, sex, BMI have slightly fewer proteins included in the testing in the UKB paper
      mutate(n_module = map(test_data, \(td) sum(td$module == i)),
             n_other = map(test_data, \(td) sum(td$module != i)),
             pos_prop_module = map(test_data, \(td) filter(td, module == i) %>%
                                     summarise(prop = sum(logp > 4.769551) / n()) %>% pull(prop)),
             pos_prop_other = map(test_data, \(td) filter(td, module != i) %>%
                                    summarise(prop = sum(logp > 4.769551) / n()) %>% pull(prop))) %>%
      mutate(ftest = map(test_data, \(td) {
        fisher.test(table(td$logp > 4.769551, td$module))
      })) %>%
      mutate(p_value = map(ftest, "p.value"),
             lower_ci = map(ftest, \(res) res[["conf.int"]][1]),
             upper_ci = map(ftest, \(res) res[["conf.int"]][2])) %>%
      select(-c(data, test_data, ftest)) %>%
      unnest(c(n_module, n_other, pos_prop_module, pos_prop_other, p_value, lower_ci, upper_ci)) %>%
      mutate(module = i)
    
  }) %>% bind_rows() %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    relocate(fdr, .after = p_value)
  
  # Compute proportions for plotting differences
  ukb_prop <- map(wgcna_ukb_assoc, \(x) {
    x %>%
      filter(logp == max(logp), .by = c(assay, trait)) %>%
      distinct(assay, trait, trait_num, logp) %>%
      count(trait_num, signif = logp > 4.769551) %>%
      mutate(prop = n / sum(n), .by = c(trait_num)) %>%
      # Fill in NAs with 0 if traits with no significant proteins at all are to be included
      # complete(trait_num, signif, fill = list(n = 0, prop = 0)) %>%
      filter(signif)
  })
  # Compute differences in proportions and add p- and FDR values
  ukb_comp <- imap(ukb_prop, \(x, i) {
    left_join(x, ukb_prop$grey %>% select(trait_num, prop) %>%
                rename(prop_grey = prop), by = "trait_num") %>%
      mutate(grey_comp = prop - prop_grey) %>%
      mutate(module = i) %>%
      left_join(tests %>% select(trait_num, p_value, fdr, module),
                by = c("trait_num", "module"))
  })

  # Make plot
  plt_dat <- map(ukb_comp, \(x) {
    x %>% select(trait_num, grey_comp, module, p_value, fdr)
  }) %>% bind_rows() %>%
    filter(module != "grey") %>%
    # Include name of trait as well
    left_join(wgcna_ukb_assoc$turquoise %>% select(trait, trait_num), by = "trait_num",
              relationship = "many-to-many")

  ukb_comp_plot <-  plt_dat %>%
    # Order traits by number
    ggplot(aes(x = module, y = reorder(trait, -trait_num))) +
    geom_tile(aes(fill = grey_comp), colour = "white") +
    # geom_text(aes(label = p_label)) +
    # geom_point(data = subset(plt_dat, fdr < 0.05), shape = "*", size = 5) +
    geom_text(aes(label = signif(fdr, 2)), data = subset(plt_dat, fdr < 0.05), size = 4, angle = 45) +
    scale_fill_gradient2(low = "skyblue2", mid = "white", high = "sienna2") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    labs(x = "Module", y = NULL, fill = "Difference in\nproportion") +
    coord_fixed() +
    theme_classic(14) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "left",
          legend.key.width = unit(0.07, "npc"))

  return(list("test_results" = tests, "plot" = ukb_comp_plot))
}

# Cluster samples based on the WGCNA MEs
run_clustering_me <- function(prot_wgcna, linearise_data = F, scale_data = F, rseed = 123, clust_distance = "manhattan", clust_method = "ward.D2", gap_b = 100, nrep = 50, verbose = F) {
  # gap_b, number of times to repeat the clustering in the gap analysis
  # nrep, number of repetitions in the stability assessment
  
  clust_in <- prot_wgcna$wgcna$MEs %>%
    # Exclude the grey module containing the non-correlating proteins
    select(-MEgrey) %>%
    # Undo log2 if desired
    (\(x) {if (linearise_data) 2**x else x})() %>%
    # Scale data if desired
    (\(x) {if (scale_data) scale(x) else x})()
  
  if (is.numeric(rseed)) {
    set.seed(rseed)
  }
  
  k <- get_optimal_k(clust_in, k_max = 15, b = gap_b, verbose = verbose,
                     distance = "manhattan", method = "ward.D2")$optimal_k
  
  stab <- clust_in %>%
    dist(method = clust_distance) %>%
    clusterboot(B = nrep, distances = T, clustermethod = hclustCBI, k = k, method = clust_method, count = verbose)
  
  # Extract the original clustering
  clust <- cutree(stab$result$result, k)
  
  clust_df <- tibble("unique_sample_name" = rownames(clust_in),
                     "cluster" = clust[unique_sample_name])
  stab_df <- count(clust_df, cluster) %>%
    arrange(cluster) %>%
    mutate(mean_ji = stab$bootmean)
  
  # Merge clusters with a mean jaccard index stability below 0.5 or a sample size below 10 to reduce the number of noise clusters
  noise_clust <- stab_df %>% filter(mean_ji < 0.5 | n < 10) %>% pull(cluster)
  # Rename other clusters in order of size
  new_clust <- stab_df %>% mutate(cluster = case_when(cluster %in% noise_clust ~ 0, T ~ cluster)) %>%
    filter(cluster != 0) %>% arrange(desc(n)) %>% mutate(cluster = 1:nrow(.))
  
  # Keep key of old vs new names to know stability of clusters
  stab_df <- stab_df %>% rename(cluster_og = cluster) %>%
    left_join(new_clust, by = c("n", "mean_ji")) %>%
    mutate(cluster = case_when(is.na(cluster) ~ 0, T ~ cluster)) %>%
    arrange(n) %>% relocate(cluster, cluster_og)
  
  # Update cluster names
  clust_df <- clust_df %>% rename(cluster_og = cluster) %>%
    left_join(stab_df %>% select(cluster_og, cluster), by = "cluster_og") %>%
    select(-cluster_og)
  
  return(list("clust_dat" = clust_df,
              "stab_dat" = stab_df))
}

# Make PCA biplot of protein modules from WGCNA
make_prot_pca_biplot <- function(prot_wgcna, prot_clust) {
  # Plot data with cluster assignments
  plt_dat <- prot_wgcna$wgcna$MEs %>%
    as.data.frame() %>%
    select(-MEgrey) %>%
    rownames_to_column("unique_sample_name") %>%
    left_join(prot_clust$clust_dat, by = "unique_sample_name") %>%
    mutate(cluster = as.character(cluster))
  
  colrs <- c("grey", brewer.pal(length(unique(plt_dat$cluster)) - 1, "Dark2")) %>%
    setNames(sort(unique(plt_dat$cluster)))
  
  plt_out <- plot_pca(plt_dat %>% select(-unique_sample_name, -cluster), pcs = 1:2, biplot = T,
                      col_vec = pull(plt_dat, cluster), scale = T, center = T,
                      col_scale = scale_colour_manual(values = colrs),
                      fill_scale = scale_fill_manual(values = colrs)) +
    labs(colour = "Proteotype", fill = "Proteotype")
  
  return(plt_out)
}

# Compare cluster memberships between serology and proteomics
comp_clust_sero_prot <- function(seroclust, protclust_dat) {
  clust_dat <- seroclust$serocluster %>% select(unique_sample_name, serocluster) %>%
    left_join(protclust_dat$clust_dat %>% rename(proteotype = cluster), by = "unique_sample_name")
  plt_dat <- clust_dat %>%
    count(serocluster, proteotype) %>%
    filter(!is.na(proteotype)) %>%
    mutate(percent = n / sum(n) * 100, .by = serocluster,
           lab = paste0(n, "\n(", signif(percent, 3), "%)"))
  
  # Plot all clusters
  plt_all <- ggplot(plt_dat, aes(x = serocluster, y = proteotype)) +
    # geom_point(aes(size = n, fill = n), shape = 21) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = lab)) +
    scale_fill_gradient(low = "white", high = "sienna2") +
    theme_classic()
  
  # Exclude cluster 0 of the proteotypes
  plt_subset <- plt_dat %>% filter(proteotype != 0) %>%
    ggplot(aes(x = serocluster, y = proteotype)) +
    # geom_point(aes(size = n, fill = n), shape = 21) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = lab)) +
    scale_fill_gradient(low = "white", high = "sienna2") +
    theme_classic()
  
  # Barplot using `comp_groups()` with seroclusters on x-axis
  n_proteotypes <- length(unique(clust_dat$proteotype)) - 1
  plt_bar <- clust_dat %>%
    mutate(across(everything(), as.character)) %$%
    comp_groups(serocluster, proteotype, p_stars = T, p_step = 0.2,
                # Change colours so that 0 is grey, matching the colours in the effects plot
                discr_palette = scale_fill_manual(values = c(
                  "0" = "grey85", brewer.pal(n_proteotypes, "Set2") %>%
                    setNames(seq(1, n_proteotypes))
                )),
                x_name = "Serocluster", y_name = "Proteotype", x_sizes = F) %>%
    pluck("plt")
  
  return(list("plt_all" = plt_all, "plt_subset" = plt_subset, "plt_bar" = plt_bar))
}

# Characterise WGCNA-based proteoclusters
characterise_proteoclust <- function(prot_wgcna, wgcna_clust, prot_dat, seroclust, crp_dat, ifn_seropos) {
  
  # Heatmap of samples in their clusters vs MEs (without the grey cluster with non-correlating proteins)
  hm_dat <- prot_wgcna$wgcna$MEs %>% select(-MEgrey) %>% as.matrix()
  # Splitting rows
  rsplit <- wgcna_clust$clust_dat$cluster[match(rownames(prot_wgcna$wgcna$MEs),
                                                wgcna_clust$clust_dat$unique_sample_name)]
  col_fun <- colorRamp2(breaks = c(min(hm_dat), 0, max(hm_dat)), colors = c("blue", "white", "red"))
  hm_out <- Heatmap(hm_dat, col = col_fun, row_split = rsplit, show_row_names = F, show_column_names = T,
                    clustering_method_rows = "ward.D2", clustering_distance_rows = "manhattan",
                    name = "ME value")
  
  # IFN classes, positive for any
  ifn_any_pos <- ifn_seropos %>% filter(str_detect(name, "IFN")) %>%
    summarise(ifn_any_pos = any(seroclass == "Positive"), .by = unique_sample_name) %>%
    mutate(ifn_any_pos = case_when(ifn_any_pos ~ "Positive", T ~ "Negative"))
  
  comp_in <- wgcna_clust$clust_dat %>%
    mutate(cluster = factor(cluster, levels = as.character(sort(unique(cluster))))) %>%
    # Add questionnaire info
    left_join(prot_dat$sinfo %>% select(unique_sample_name, age_group, sex,
                                        infected, vaccinated, symptoms, inf_diff, vacc_diff),
              by = "unique_sample_name") %>%
    mutate(infected = as.character(infected), vaccinated = as.character(vaccinated),
           symptoms = factor(symptoms, levels = c("None", "Smell", "Taste", "Smell+Taste"))) %>%
    # Add serocluster memberships
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster),
              by = "unique_sample_name") %>%
    # Add CRP classification
    left_join(crp_dat %>% select(unique_sample_name, crp_class = class),
              by = "unique_sample_name") %>%
    # Add IFN AAb classes (positive for any anti-IFN AAb)
    left_join(ifn_any_pos, by = "unique_sample_name")
  
  # Pairwise comparisons of clusters
  comp_out <- sapply(colnames(comp_in)[-c(1:2)], \(x) {
    comp_groups(as.formula(paste0(x, " ~ cluster")), comp_in,
                discr_palette = case_when(x == "serocluster" ~ "Purples",
                                          x == "crp_class" ~ "Greens",
                                          T ~ "Set2"),
                p_stars = T, p_step = 0.2, x_name = "Proteotype", x_sizes = F,
                y_name = case_when(x == "serocluster" ~ "Serocluster",
                                   x == "crp_class" ~ "CRP level",
                                   x == "symptoms" ~ "Symptoms",
                                   x == "ifn_any_pos" ~ "Positive for any IFN",
                                   T ~ x))
  }, simplify = F, USE.NAMES = T)
  
  comp_p <- map(comp_out, "p.value")
  comp_plt <- map(comp_out, "plt")
  
  # Comparisons in the opposite direction (trait on x axis)
  comp_out_rev <- sapply(comp_in %>% select(-c(unique_sample_name, cluster, inf_diff, vacc_diff)) %>%
                           colnames(), \(x) {
    comp_groups(as.formula(paste0("cluster ~ ", x)), comp_in %>% filter(!is.na(vaccinated)),
                discr_palette = scale_fill_manual(values = c("grey85", brewer.pal(length(unique(comp_in$cluster)) - 1, "Set2")) %>% setNames(levels(comp_in$cluster))),
                p_stars = F, p_step = 0.2, x_sizes = T, y_name = "Proteotype",
                x_name = case_when(x == "serocluster" ~ "Serocluster",
                                   x == "crp_class" ~ "CRP level",
                                   x == "vaccinated" ~ "Vaccine doses"))
  }, simplify = F, USE.NAMES = T)
  comp_p_rev <- map(comp_out_rev, "p.value")
  comp_plt_rev <- map(comp_out_rev, "plt")
  
  # Table with overall comparisons
  # Kruskal-Wallis test for numeric variables and Fisher exact test for categorical variables
  pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
      # For numeric variables, perform a standard 2-sample t-test
      p <- kruskal.test(y ~ g)$p.value
    } else {
      # For categorical variables, perform a chi-squared test of independence
      p <- fisher.test(table(y, g), simulate.p.value = T)$p.value
    }
    # Format the p-value, using scientific notation for small values
    # The initial empty string places the output on the line below the variable label.
    c("", if_else(p < 0.01,
                  as.character(formatC(p, format = "e", digits = 2)),
                  as.character(signif(p, digits = 2))))
  }
  tbl1 <- table1(~ age_group + sex + infected + vaccinated +
                   symptoms + serocluster + crp_class + ifn_any_pos | cluster,
                 data = comp_in, overall = F, extra.col = list(`P-value` = pvalue))
  
  return(list("heatmap" = hm_out, "p_value" = comp_p, "plot" = comp_plt,
              "p_value_rev" = comp_p_rev, "plot_rev" = comp_plt_rev, "table" = tbl1))
}

# Calculate effect sizes (differences in mean) of proteome clusters/phenotypes (cluster vs reference or cluster vs cluster)
compute_clust_effect <- function(x_in, prot_clust, prot_wgcna, clust_vs_rest = T) {
  # x, data frame containing the columns unique_sample_name (for binding together with cluster assignments), a column named 'name' with the name of the analyte (protein, module, whatever you want to compare), and a column named value with the values for the analytes
  
  # Compute 95% confidence intervals
  d_ci <- x_in %>%
    left_join(prot_clust$clust_dat, by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, cluster)) %>%
    nest(.by = name) %>%
    mutate(diff_ci = map(data, \(x) {
      
      comp_out <- if (clust_vs_rest) {
        # Cluster vs rest
        map(x %>% filter(cluster != 0) %>% pull(cluster) %>% unique(), \(y) {
          x1 <- x %>% filter(cluster == y)
          x2 <- x %>% filter(cluster != y)
          
          tt <- t.test(x1$value, x2$value)
          mean_diff <- mean(x1$value) - mean(x2$value)
          
          data.frame(comparison = y,
                     mean_diff = mean_diff,
                     lower = tt$conf.int[1],
                     upper = tt$conf.int[2])
        })
        
      } else {
        # Cluster vs cluster
        map(combn(sort(unique(prot_clust$clust_dat$cluster)), 2, simplify = F), \(y) {
          x1 <- x %>% filter(cluster == y[1])
          x2 <- x %>% filter(cluster == y[2])
          
          tt <- t.test(x1$value, x2$value)
          mean_diff <- mean(x1$value) - mean(x2$value)
          
          data.frame(comparison = paste(y, collapse = "-"),
                     mean_diff = mean_diff,
                     lower = tt$conf.int[1],
                     upper = tt$conf.int[2])
        })
      }
      
      return(bind_rows(comp_out))
    })) %>%
    select(name, diff_ci) %>%
    unnest(diff_ci) %>%
    mutate(module = prot_wgcna$wgcna$colors[name])
  
  return(d_ci)
}

# Visualise protein differences between clusters
make_prot_clust_effect_plots <- function(plt_dat, aes_pointrange = aes(x = mean_diff, y = name, xmin = lower, xmax = upper, colour = comparison, group = desc(comparison)),
                                         spread = 0.5, point_size = 1, linewidth = 1,
                                         y_min_grid = seq(1.5, length(unique(plt_dat$name)) - 0.5),
                                         min_grid_lty = "solid", min_grid_lwd = 0.5, min_grid_col = "grey",
                                         x_name = "Effect size", y_name = NULL, col_name = "Proteotype") {
  # spread, number to set as the `width` argument to position_dodge (spread of cluster lines within each protein)
  # y_min_grid, y axis minor grid drawn by `geom_hline()`. By default the lines are drawn between each element on the y axis i.e. at 1.5, 2.5, 3.5, ..., but the breaks can be adjusted by providing a nimeric vector with the breaks. If no grid is desired, set to NULL.
  
  # Plot proteins
  ci_plot <- plt_dat %>%
    ggplot() +
    geom_hline(yintercept = y_min_grid, colour = min_grid_col, linetype = min_grid_lty,
               linewidth = min_grid_lwd) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_pointrange(aes_pointrange, size = point_size, position = position_dodge(width = spread),
                    linewidth = linewidth) +
    scale_colour_brewer(palette = "Set2") +
    scale_x_continuous(limits = center_limits()) +
    labs(x = x_name, y = NULL, colour = col_name) +
    theme_classic(14)
  
  return(ci_plot)
}

save_proteotype_fig <- function(prot_me_pca, prot_clust_char, prot_me_effect_plts,
                                prot_clust_effect, prot_seroclust_plts, wgcna_prot_overlap) {
  filename <- make_fn("protein_clustering", "prot_clustering_fig.pdf")
  
  # Add example proteins for each module as a blank plot with labels
  prot_ex_df <- prot_clust_effect %>%
    # Count proteins with significant difference in mean (based on 95% CI)
    add_count(name, n_signif = !(lower < 0 & upper > 0)) %>%
    # The .drop argument of add_count() does nothing and is deprecated, meaning that if any protein has no
    # significant comparisons or only significant comparisons there is no 0 added for the class that has
    # zero counts. Add those manually for those with zero significant comparisons to not lose when filtering
    mutate(n = case_when(!n_signif & n == 5 ~ 0, T ~ n)) %>%
    mutate(n_signif = case_when(n == 0 ~ T, T ~ n_signif)) %>%
    filter(n_signif) %>%
    distinct(name, .keep_all = T) %>%
    mutate(name_short = str_remove(name, "_\\d$")) %>%
    # Add variable indicating if the protein is in one of the protein vs serocluster boxplots (if so prioritise)
    # Also prioritise some selected proteins (CC16 (SCGB1A1) and MMP10)
    # Finally prioritise proteins with two measurements where both end up in the same WGCNA module
    # If in none of the categories, sort by number of significant diferences and then magnitude of difference
    mutate(in_boxplt = as.numeric(name_short %in% (prot_seroclust_plts$wilcox_box$protein %>%
                                                     str_remove("_\\d$"))),
           prio_prot = as.numeric(name_short %in% c("SCGB1A1", "MMP10")),
           overlap = as.numeric(name_short %in% (wgcna_prot_overlap$data %>%
                                                   filter(same_module_across_technologies) %>%
                                                   pull(assay)))) %>%
    # Keep only one of the two if multiple
    distinct(name_short, .keep_all = T) %>%
    slice_max(order_by = tibble(prio_prot, in_boxplt, overlap, n, abs(mean_diff)), n = 3, by = module) %>%
    # Compute y positions based on module positions
    mutate(module = str_to_sentence(module),
           module = factor(module, levels = rev(append(setdiff(sort(unique(module)), "Grey"), "Grey"))),
           module_text = module,
           module = factor(module, levels = sort(unique(module))) %>%
             as.integer(),
           module = module + rep(c(-0.25, 0, 0.25), length(unique(module))),
           x = 0) %>%
    select(name_short, module, module_text, x, in_boxplt, prio_prot, overlap, n, mean_diff)
  
  # Make a separate plot with the labels to not mess with the x-axes in the effect size plot
  prot_ex_plt <- ggplot(prot_ex_df) +
    # Make a point plot with discrete axis to get the same axis as the effect size plot
    geom_point(aes(x = x, y = module_text), alpha = 0) +
    geom_text(aes(x = x, y = module, label = name_short)) +
    labs(title = "Example\nproteins") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  
  effect_plt <- wrap_plots(prot_me_effect_plts +
                             theme(axis.text.y = element_text(size = 14)),
                           prot_ex_plt, widths = c(3, 1))
  
  # Bar plots
  bar_plt <- wrap_plots(prot_clust_char$plot_rev$serocluster,
                        prot_clust_char$plot_rev$crp_class,
                        prot_clust_char$plot_rev$vaccinated, nrow = 1) +
    plot_layout(guides = "collect", axis_titles = "collect_y") &
    theme_classic(13)
  
  fig <- list(prot_me_pca, effect_plt, bar_plt) %>%
    # wrap_plots(nrow = 1, widths = c(4, 1)) &
    wrap_plots(design = paste0(paste0(rep("AAAA\n", 4), collapse = ""),
                               paste0(rep("BBBB\n", 4), collapse = ""),
                               paste0(rep("CCCC\n", 2), collapse = ""))) &
    coord_cartesian(clip = "off") &
    theme(legend.title.position = "left",
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0))
  ggsave(filename = filename, plot = fig, height = 12, width = 8)
  
  return(filename)
}
