# Interferon target functions

# Data with just interferons for convenience
make_ifn_data <- function(serol_in) {
  dat <- serol_in
  dat$sinfo <- dat$sinfo %>% filter(class == "Sample")
  dat$mfi <- dat$mfi %>%
    filter(unique_sample_name %in% dat$sinfo$unique_sample_name) %>%
    # Include data for anti-human IgG and anti-human IgM as well for a baseline of Ab levels
    select(unique_sample_name, starts_with("IFN"), anti_human_igg, anti_human_igm)
  dat$sinfo <-  dat$sinfo %>% filter(unique_sample_name %in% dat$mfi$unique_sample_name)
  
  dat$binfo <- dat$binfo %>%
    filter(unique_antigen_name %in% colnames(dat$mfi)) %>%
    # Generate the long form protein names from the acronyms
    mutate(protein_name = str_replace_all(unique_antigen_name, c(
      "IFN" = "Interferon ", "A" = "alpha ", "B" = "beta ",
      "G" = "gamma ", "L" = "lambda ", "W" = "omega ", "R" = "receptor "
    )))
  
  return(dat[c("mfi", "sinfo", "binfo", "comment")])
}

# Classify seropositivity of interferons
classify_ifn <- function(ifn_in, n_sd = 12) {
  # Parameter for cutoff
  neg_prop <- 1 - sum(ifn_in$sinfo$vaccinated > 0 |
                      ifn_in$sinfo$infected > 0, na.rm = T) /
    nrow(ifn_in$sinfo)
  
  seropos_class <- ifn_in$mfi %>%
    pivot_longer(cols = -unique_sample_name) %>%
    mutate(seroclass = density_cutoff(value, n_sd, neg_prop)$classes,
           cutoff = density_cutoff(value, n_sd, neg_prop)$cutoff,
           seroclass = case_when(seroclass ~ "Positive", T ~ "Negative"),
           .by = name) %>%
    select(-value)
    
  return(seropos_class)
}

# Plot interferon levels in each cluster
make_ifn_clust_plots <- function(ifn_in, ifn_seropos, seroclust) {
  
  plt_dat <- ifn_in$mfi %>%
    # Add cluster membership
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster),
              by = "unique_sample_name") %>%
    filter(!is.na(serocluster)) %>%
    pivot_longer(cols = -c(unique_sample_name, serocluster)) %>%
    # Add seropositivity
    left_join(ifn_seropos, by = c("unique_sample_name", "name")) %>%
    # Reorder plots
    mutate(name = case_when(str_detect(name, "igg") ~ "Anti-human IgG",
                            str_detect(name, "igm") ~ "Anti-human IgM",
                            T ~ name)) %>%
    mutate(name = factor(name, levels = c(
      sort(ifn_in$binfo$unique_antigen_name), "Anti-human IgG", "Anti-human IgM"
    )))
  
  boxplots <- ggplot(plt_dat, aes(x = serocluster, y = value)) +
    geom_hline(aes(yintercept = cutoff), colour = "grey50", linetype = "dotted") +
    # Beeswarm points + boxplot
    ggbeeswarm::geom_beeswarm(aes(colour = seroclass), alpha = 0.5) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    scale_colour_manual(values = c("Negative" = "grey80", "Positive" = "blue")) +
    facet_wrap2(~ name, nrow = 3) +
    labs(x = "Serocluster", y = "AAb levels [AU]", colour = "Seropositivity") +
    theme_classic()
  
  # Bar plots with pairwise Fisher test
  # Use the comp_groups function to perform the tests
  pvalues <- plt_dat %>%
    select(name, serocluster, seroclass) %>%
    group_by(name) %>% nest() %>%
    mutate(p = map(data, \(x) {
      # NA if only one level
      if (length(unique(x$seroclass)) == 1) return(data.frame(clust2 = NA, clust1 = NA, pval = NA, name = name))
      
      comp_groups(seroclass ~ serocluster, x, p_adj_method = "none")$p.value %>%
        as.data.frame() %>%
        rownames_to_column("clust2") %>%
        pivot_longer(cols = -clust2, names_to = "clust1", values_to = "pval") %>%
        filter(!is.na(pval)) %>%
        mutate(name = name)
    })) %>%
    pull(p) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    relocate(name, clust1, clust2)
  
  # Make another data frame with plotting data, use for computing p-value positions
  plt_dat2 <- plt_dat %>%
    # Turn cluster names into a factor to include 0 in the counts
    mutate(serocluster = factor(serocluster, levels = as.character(1:4))) %>%
    count(name, serocluster, seroclass, .drop = F) %>%
    mutate(percent = n / sum(n) * 100, .by = c(name, serocluster)) %>%
    filter(seroclass == "Positive")
  
  # Pick out FDR values below 0.05 for plot, calculate y positions of labels
  plabels <- pvalues %>%
    filter(fdr < 0.05) %>%
    mutate(y = max(plt_dat2$percent) +                                # Start from highest value
             0.1 * diff(range(plt_dat2$percent)) +                    # Add 10% of range to be above max
             0.2 * diff(range(plt_dat2$percent)) * (seq(1, n()) - 1), # Add 20% per extra annotation
           .by = name) %>%
    mutate(label = paste0("italic('p') == ", signif(fdr, 3)))
  
  barplots <- plt_dat2 %>%
    ggplot() +
    geom_bar(aes(x = serocluster, y = percent), stat = "identity") +
    geom_signif(aes(xmin = clust1, xmax = clust2, y_position = y, annotations = label),
                data = plabels, manual = T, parse = T) +
    facet_wrap2(~ name, nrow = 3) +
    labs(x = "Serocluster", y = "Frequency in cluster [%]") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_classic()
  
  return(list("boxplots" = boxplots, "barplots" = barplots, "pvalues" = pvalues))
}

# Make table of interferon positive percentages in total and in clusters
make_ifn_table <- function(ifn_in, ifn_seropos, seroclust) {
  
  # Overall table
  tbl_overall <- ifn_seropos %>%
    summarise(percent_pos = sum(seroclass == "Positive") / n() * 100,
              # Add confidence intervals
              ci_pos = paste(round(prop_ci(percent_pos / 100, n()) * 100, 1), collapse = ","),
              "All [%]" = paste0(round(percent_pos, 1), " (", ci_pos, ")"),
              .by = name) %>%
    select(name, percent_pos, "All [%]")
  
  # Cluster table
  tbl_cluster <- ifn_seropos %>%
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster), by = "unique_sample_name") %>%
    summarise(n_pos = sum(seroclass == "Positive"),
              percent_pos_cluster = n_pos / n() * 100,
              ci_pos_cluster = paste(round(prop_ci(percent_pos_cluster / 100, n()) * 100, 1), collapse = ","),
              
              # percent_pos_total = n_pos / nrow(ifn_in$mfi) * 100,
              .by = c(name, serocluster)) %>%
    mutate(pos = paste0(round(percent_pos_cluster, 1), " (", ci_pos_cluster, ")")) %>%
    arrange(name, serocluster) %>%
    pivot_wider(id_cols = "name", names_from = "serocluster", values_from = "pos")
  
  tbl_combined <- tbl_overall %>%
    left_join(tbl_cluster, by = "name") %>%
    left_join(ifn_in$binfo %>% select(unique_antigen_name, protein_name), by = c("name" = "unique_antigen_name")) %>%
    # Look at only interferons
    filter(str_detect(name, "IFN")) %>%
    rename("Gene" = "name", "Protein" = "protein_name") %>%
    relocate("Protein", .after = 1) %>%
    arrange(desc(percent_pos)) %>%
    select(-percent_pos)
  
  return(tbl_combined)
}

# Analyse age, sex and SARS-CoV-2 serostatus in relation to IFN levels
check_ifn_assoc <- function(ifn_in, ifn_seropos, cov_seropos) {
  # Bind together age, sex and symptom information with seropositivity information (use MFI to get Abs in there too)
  d <- ifn_in$mfi %>%
    left_join(ifn_in$sinfo %>% select(unique_sample_name, age_group, sex, symptoms) %>%
                mutate(symptoms = case_when(symptoms != "None" ~ "Yes", T ~ "No")),
              by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, age_group, sex, symptoms)) %>%
    # IFN seropositivity
    left_join(ifn_seropos %>% select(-cutoff), by = c("unique_sample_name", "name")) %>%
    # Anti-SARS-CoV-2 seropositivity, convert TRUE/FALSE to "Positive"/"Negative"
    left_join(cov_seropos$per_cov_prot, by = "unique_sample_name") %>%
    mutate(across(where(is.logical), \(x) {case_when(x ~ "Positive", T ~ "Negative")})) %>%
    # Reorder plots
    mutate(name = case_when(str_detect(name, "igg") ~ "Anti-human IgG",
                            str_detect(name, "igm") ~ "Anti-human IgM",
                            T ~ name)) %>%
    mutate(name = factor(name, levels = c(
      sort(ifn_in$binfo$unique_antigen_name), "Anti-human IgG", "Anti-human IgM"
    )))
  
  # Little function for the plots
  bar_pval_plot <- function(dat_in, comp_var) {
    plt_dat <- dat_in %>%
      group_by(grp = get(comp_var), name, seroclass) %>%
      filter(!is.na(grp)) %>%
      summarise(n = n(), .groups = "keep") %>%
      ungroup(seroclass) %>%
      mutate(percent = n / sum(n) * 100) %>%
      ungroup()
    
    # Compute p-values with fisher test for the grouping variable
    grp_pval <- dat_in %>%
      summarise(p = tryCatch(table(get(comp_var), seroclass) %>% fisher.test() %>% pluck("p.value"),
                             # If not enough observations an error is thrown. If so, return NA as no test is run
                             error = \(err) {NA}),
                .by = name) %>%
      mutate(fdr = p.adjust(p, method = "fdr")) %>%
      mutate(label = case_when(!is.na(fdr) ~ fdr %>% signif(3) %>% paste0("italic('p') == ", .),
                               T ~ ""))
    
    # Pre-filter to use max and range of positive values for annotations
    plt_dat <- plt_dat %>% filter(seroclass == "Positive")
    
    plt <- plt_dat %>%
      ggplot(aes(x = grp, y = percent)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = label),
                x = mean(as.numeric(as.factor(unique(plt_dat$grp)))),
                y = max(plt_dat$percent) + 0.1 * diff(range(plt_dat$percent)),
                parse = T, data = grp_pval) +
      facet_wrap(~ name, nrow = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(x = comp_var %>% str_replace_all("_", " ") %>% str_to_sentence(), y = "Percent") +
      theme_classic()
    
    return(list("plot" = plt, "pvalue" = grp_pval, "dat" = dat_in))
  }
  
  # Make plots
  plots_n_pvals <- sapply(c("age_group", "sex", "symptoms", "N+", "S+"), \(x) {bar_pval_plot(d, x)},
                          simplify = F, USE.NAMES = T)
  
  return(plots_n_pvals)
}

make_ifn_assoc_summary <- function(ifn_assoc) {
  # Extract all p-values
  pvals <- imap(ifn_assoc, \(data, index) {
    data$pvalue %>% mutate(var = index)
  }) %>%
    bind_rows() %>%
    select(-label) %>%
    # Rename and reorder variable names for plot
    mutate(var = str_replace(var, "age_group", "age") %>% str_to_sentence()) %>%
    mutate(var = factor(var, levels = c("Age", "Sex", "Symptoms", "N+", "S+"))) %>%
    mutate(is_ifna = case_when(str_detect(name, "IFNA\\d") ~ "Yes", T ~ "No"))
  
  # p-value plot
  p_plot <- pvals %>%
    filter(str_detect(name, "IgG|IgM", negate = T)) %>%
    mutate(is_ifna = case_when(is_ifna == "Yes" ~ "IFNA", is_ifna == "No" ~ "Other")) %>%
    ggplot(aes(x = var, y = -log10(fdr), colour = is_ifna)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
    geom_beeswarm(size = 3, alpha = 0.75, cex = 1.5) +
    scale_colour_brewer(palette = "Paired", direction = -1) +
    labs(x = NULL, y = "-log10(P-value)", colour = NULL) +
    theme_classic(16) +
    theme(legend.position = "inside", legend.position.inside = c(0.2, 0.85))
  
  # p-value table
  p_table <- pvals %>%
    pivot_wider(id_cols = name, names_from = var, values_from = p)
  
  return(list("plot" = p_plot, "table" = p_table))
}

# Interferon repsonse at different time points 
plot_ifn_response <- function(ifn_dat, serol_seropos, ifn_seropos) {
  
  # Number of positive individuals over time
  plts <- ifn_seropos %>%
    filter(str_detect(name, "IFN")) %>%
    select(-cutoff) %>%
    left_join(ifn_dat$sinfo %>% select(unique_sample_name, inf_diff), by = "unique_sample_name") %>%
    rename(months = inf_diff) %>%
    # Look at seropositive individuals, colour by N seroclass
    filter(!is.na(months), seroclass == "Positive") %>%
    left_join(serol_seropos$per_cov_prot %>% select(unique_sample_name, `N+`), by = "unique_sample_name") %>%
    count(name, seroclass, months, `N+`) %>%
    ggplot(aes(x = months, y = n, fill = `N+`)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ name) +
    theme_classic()
  
  return(plts)
}

# Interferon associations with proteins
check_ifn_prot_assoc <- function(ifn_seropos, prot_dat) {
  # Logistic regression of seroclass vs protein measurements
  test_results <- ifn_seropos %>%
    rename(ifn = name) %>%
    filter(str_detect(ifn, "IFN")) %>%
    select(unique_sample_name, ifn, seroclass) %>%
    left_join(prot_dat$npm, by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, ifn, seroclass),
                 names_to = "protein", values_to = "npm") %>%
    filter(!is.na(npm)) %>%
    mutate(seroclass = as.factor(seroclass)) %>%
    group_by(ifn, protein) %>% nest() %>%
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
    # unnest(p_value) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr"), .by = ifn)
  
  return(test_results)
}

# Make heatmap of IFN protein associations
make_ifn_prot_hms <- function(ifn_prot_assoc, prot_ab_assoc) {
  # Plotting data
  plt_dat <- ifn_prot_assoc %>%
    mutate(fill_p = case_when(p_value < 0.05 & estimate > 0 ~ "Positive",
                              p_value < 0.05 & estimate < 0 ~ "Negative",
                              p_value > 0.05 ~ "None"),
           fill_fdr = case_when(fdr < 0.05 & estimate > 0 ~ "Positive",
                                fdr < 0.05 & estimate < 0 ~ "Negative",
                                fdr > 0.05 ~ "None")) %>%
    select(protein, ifn, fill_p, fill_fdr, p_value, fdr) %>%
    # Reorder IFNs so the numbers are in numeric order and not alphabetic
    mutate(ifn_char = str_extract(ifn, "[:alpha:]+"),
           ifn_num = str_extract(ifn, "\\d+$") %>% as.numeric()) %>%
    arrange(ifn_char, ifn_num) %>%
    mutate(ifn = factor(ifn, levels = unique(ifn)))
  
  # Keep only proteins measured in both assays that have some common IFN hits in the nominal p heatmap
  # Also for IFNs, keep only if any protein has a hit in both methods
  any_both <- plt_dat %>%
    filter(str_detect(protein, "_\\d$")) %>%
    separate_wider_delim(protein, delim = "_", names = c("protein", "num")) %>%
    summarise(any_both = (any(fill_p[num == 1] == "Positive") & any(fill_p[num == 2] == "Positive")) |
                (any(fill_p[num == 1] == "Negative") & any(fill_p[num == 2] == "Negative")),
              .by = c(protein, ifn)) %>%
    filter(any_both)
  
  # Add some proteins from the protein Ab association to make it a protein vs circulating Abs plot
  any_both <- add_row(any_both, protein = c("CD200R1", "CXCL1", "CXCL8", "CXCL9", "TREM2"))
  
  # Prepare data for nominal p plot
  test_hm_p <- plt_dat %>%
    filter(protein %in% paste0(rep(unique(any_both$protein), each = 2), c("_1", "_2")),
           ifn %in% unique(any_both$ifn)) %>%
    # Reorder proteins to start at the top for protein splits
    arrange(protein) %>%
    mutate(protein = factor(protein, levels = unique(protein))) %>%
    # Keep IFN name without number for facetting
    select(protein, ifn, fill_p, ifn_char) %>%
    # Get the values from the prot ab assoc data
    bind_rows(
      prot_ab_assoc$heatmap_p$data %>%
        filter(protein %in% paste0(rep(unique(any_both$protein), each = 2), c("_1", "_2"))) %>%
        rename(ifn = ab) %>%
        mutate(ifn_char = case_when(str_detect(ifn, "Na|Nc") ~ "Anti-N",
                                    str_detect(ifn, "RBD|S1|S1S2") ~ "Anti-S")) %>%
        select(protein, ifn, fill_p, ifn_char)
    ) %>%
    # Make a variable to facet by each protein
    mutate(prot_split = str_remove(protein, "_\\d$"))
  
  # Use the ComplexHeatmap package for easier splitting by protein and IFN
  col_map <- structure(c("sienna2", "white", "skyblue2"), names = c("Positive", "None", "Negative"))
  test_hm_p <- test_hm_p %>%
    pivot_wider(id_cols = "protein", names_from = "ifn", values_from = "fill_p") %>%
    # Arrange in proteins in descending order (factor level order) to start alphabetically from the top
    arrange(protein) %>%
    column_to_rownames("protein") %>%
    as.matrix() %>%
    Heatmap(
      column_split = test_hm_p %>% distinct(ifn, ifn_char) %>% pull(ifn_char),
      row_split = test_hm_p %>% distinct(protein, prot_split) %>% pull(prot_split),
      gap = unit(3, "pt"), column_title = NULL, show_row_names = F, row_title_rot = 0,
      col = col_map, name = "Association", rect_gp = gpar(col = "grey", lwd = 0.5),
      heatmap_legend_param = list(labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 15))
    )
  
  # Make heatmap of results
  # Split it into multiple parts for visibility
  splits <- ceiling(seq_along(unique(ifn_prot_assoc$protein)) /
                      ceiling(length(unique(ifn_prot_assoc$protein)) / 5)) %>%
    setNames(sort(unique(ifn_prot_assoc$protein)))
  
  # Keep only IFNs with any potential hits for FDR heatmap
  ifns <- ifn_prot_assoc %>%
    count(ifn, signif = fdr < 0.05) %>%
    filter(signif) %>%
    pull(ifn)
  
  test_hm_fdr <- plt_dat %>%
    filter(ifn %in% ifns) %>%
    mutate(facet = splits[protein]) %>%
    # Reorder proteins to start at the top
    arrange(facet, desc(protein)) %>%
    mutate(protein = factor(protein, levels = unique(protein))) %>%
    select(protein, ifn, fill_fdr, facet) %>%
    ggplot(aes(x = ifn, y = protein, fill = fill_fdr)) +
    geom_tile(colour = "grey", linewidth = 0.1) +
    facet_grid2(cols = vars(facet), scales = "free", drop = T, independent = "all") +
    labs(x = NULL, y = NULL, fill = "Association") +
    scale_fill_manual(values = c("Positive" = "sienna2", "Negative" = "skyblue2", "None" = "white")) +
    theme_minimal(14) +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list("heatmap_p" = test_hm_p, "heatmap_fdr" = test_hm_fdr, "proteins_p" = any_both))
}

# Figure of IFN associations
save_ifn_summary <- function(ifn_assoc_summary) {
  fn <- make_fn("interferon_analysis", "ifn_summary.pdf")
  ggsave(fn, ifn_assoc_summary$plot, width = 4.5, height = 4)
}

