# Serology target functions

# Boxplots of serology results
make_serol_boxplots <- function(serol_in) {
  # Pick out samples
  mfi <- serol_in$mfi %>% filter(unique_sample_name %in%
                                   (serol_in$sinfo %>% filter(class == "Sample") %>% pull(unique_sample_name)))
  sinfo <- serol_in$sinfo %>% filter(unique_sample_name %in% mfi$unique_sample_name)
  
  # Make boxplots
  boxplt <- mfi %>%
    left_join(sinfo %>% select(unique_sample_name, region, vacc_inf), by = "unique_sample_name") %>%
    filter(!is.na(vacc_inf)) %>%
    pivot_longer(cols = -c(unique_sample_name, vacc_inf, region)) %>%
    # Reorder and rename groups for consistency
    mutate(vacc_inf = str_remove(vacc_inf, " "),
           vacc_inf = factor(vacc_inf, levels = c("V-I-", "V+I-", "V-I+", "V+I+"))) %>%
    # One plot per Ab, self-reported vaccination/infection status on x axis, split boxes by region
    nest_by(name) %>%
    mutate(plt = list(
      ggplot(data, aes(x = vacc_inf, y = value, colour = region)) +
        geom_beeswarm(dodge.width = 0.75, size = 2) +
        geom_boxplot(aes(fill = region), width = 0.2, alpha = 0,
                     colour = "black", linewidth = 0.7, position = position_dodge(width = 0.75),
                     outlier.shape = NA, show.legend = F) +
        labs(x = "Self-reported immune category", y = name, colour = "Region") +
        scale_colour_brewer(palette = "Paired") +
        theme_classic(20)
      
    )) %>%
    select(-data)
  
  return(boxplt %>% pull(plt) %>% setNames(boxplt$name))
}

# Seropositivity classification
classify_serology <- function(serol_in, n_sd = 6) {
  # Pick out samples
  mfi <- serol_in$mfi %>% filter(unique_sample_name %in%
                                   (serol_in$sinfo %>% filter(class == "Sample") %>% pull(unique_sample_name)))
  sinfo <- serol_in$sinfo %>% filter(unique_sample_name %in% mfi$unique_sample_name)
  
  neg_prop_n <- 1 - sum(sinfo$infected > 0, na.rm = T) / nrow(sinfo)
  neg_prop_s <- 1 - sum(sinfo$infected > 0 | sinfo$vaccinated > 0, na.rm = T) / nrow(sinfo)
  
  seropos <- mfi %>% get_cov() %>%
    pivot_longer(cols = -unique_sample_name) %>%
    group_by(name) %>% nest() %>%
    mutate(seroclass = map(data, \(x) {
      cls <- density_cutoff(x$value, n_sd, ifelse(str_detect(name, regex("anti_n", ignore_case = T)), neg_prop_n, neg_prop_s))$classes
      
      return(case_when(cls ~ "Positive", T ~ "Negative"))
    })) %>%
    unnest(c(data, seroclass)) %>%
    ungroup()
  
  # Classes where antibodies targeting the same protein are bunched together
  seropos_small <- seropos %>%
    pivot_wider(id_cols = "unique_sample_name", names_from = "name", values_from = "seroclass") %>%
    # Seropositive if more than one is positive
    summarise(`N+` = case_when(sum(anti_nc == "Positive", anti_na == "Positive") > 1 ~ T, T ~ F),
              `S+` = case_when(sum(anti_s1 == "Positive", anti_s1s2 == "Positive", anti_rbd == "Positive") > 1 ~ T, T ~ F),
              .by = "unique_sample_name")
  
  return(list("per_ab" = seropos, "per_cov_prot" = seropos_small))
}

# Serology summary table
make_serol_table <- function(serol_in, seropos_in) {
  
  seropos <- seropos_in$per_ab %>% summarise(percent = signif(sum(seroclass == "Positive") / n() * 100, digits = 3),
                                             ci = paste(signif(prop_ci(percent / 100, n()) * 100, digits = 3),
                                                        collapse = ","),
                                             .by = name)
  
  tbl <- serol_in$binfo %>%
    filter(str_detect(unique_antigen_name, "anti_(n|s|rbd)")) %>%
    # Add seropositivity information
    left_join(seropos, by = c("unique_antigen_name" = "name")) %>%
    # Rename for consistent names across figures and tables
    mutate(Acronym = rename_serol(unique_antigen_name)) %>%
    select("Acronym", "SARS-CoV-2 protein" = "Antigen description",
           "Source" = "Provider", "% positive" = "percent",
           "95% CI" = "ci") %>%
    arrange(Acronym)
  
  return(tbl)
}

# Hierarchical clustering of serology data
run_serol_hclust <- function(serol_in) {
  sinfo <- serol_in$sinfo
  binfo <- serol_in$binfo
  mfi <- serol_in$mfi %>% get_cov()
  
  # Samples to include
  s_incl <- sinfo %>%
    filter(class == "Sample") %>%
    pull(unique_sample_name)
  sinfo <- sinfo[match(s_incl, sinfo$unique_sample_name), ]
  
  # Scaling of data
  mfi <- mfi[na.omit(match(s_incl, mfi$unique_sample_name)), ] %>%
    column_to_rownames("unique_sample_name") %>%
    scale(center = T, scale = T)
  
  sinfo <- sinfo %>% filter(unique_sample_name %in% rownames(mfi))
  
  # Distance and method for clustering
  clust_dist <- "manhattan"
  clust_method <- "ward.D2"
  n_clust <- 4
  
  # Perform stability test
  set.seed(123)
  bs <- clusterboot(mfi, B = 100, clustermethod = boot_clustfun, count = F,
                    dist_method = clust_dist, cl_method = clust_method, n_cl = n_clust)
  
  # Get results from first run as clustering to look at
  hcl <- bs$result$result
  
  clust <- cutree(hcl, k = n_clust)
  
  serocluster <- sinfo %>% select(unique_sample_name, vacc_inf) %>%
    mutate(og_serocluster = as.character(clust[match(unique_sample_name, names(clust))]))
  
  # Rename cluster numbers to be in size order
  serocluster_rewrite <- serocluster %>% count(og_serocluster) %>%
    arrange(desc(n)) %>% mutate(serocluster = as.character(1:nrow(.)))
  
  serocluster <- left_join(serocluster, serocluster_rewrite %>% select(og_serocluster, serocluster),
                           by = "og_serocluster")
  
  # Vaccination/Infection status associated with each cluster
  clust_vi <- serocluster %>%
    group_by(serocluster) %>%
    count(vacc_inf) %>%
    arrange(serocluster) %>%
    filter(n == max(n))
  
  serocluster <- left_join(serocluster, clust_vi %>% select(serocluster, serocluster_vi = vacc_inf),
                           by = "serocluster") %>%
    mutate(serocluster_vi = factor(serocluster_vi, levels = c("V- I-", "V+ I-", "V- I+", "V+ I+")),
           vacc_inf = factor(vacc_inf, levels = levels(serocluster_vi)))
  
  # Keep the old cluster labels together with stability
  serocluster_rewrite$mean_ji <- bs$bootmean[as.integer(serocluster_rewrite$og_serocluster)]
  
  hcl_out <- list("clustering" = hcl,
                  "serocluster" = serocluster %>% select(-og_serocluster),
                  "cluster_stab" = serocluster_rewrite,
                  "data" = mfi %>% as.data.frame() %>% rownames_to_column("unique_sample_name"))
  return(hcl_out)
}

# Summary table of clusters
make_seroclust_table <- function(seroclust, serol_in, seropos, ifn_seropos) {
  
  # IFN classes, positive for any
  ifn_any_pos <- ifn_seropos %>% filter(str_detect(name, "IFN")) %>%
    summarise(ifn_any_pos = any(seroclass == "Positive"), .by = unique_sample_name) %>%
    mutate(ifn_any_pos = case_when(ifn_any_pos ~ "Positive", T ~ "Negative"))
  
  # Table in the table1 style
  d <- seroclust$serocluster %>%
    left_join(serol_in$sinfo %>% select(unique_sample_name, age_group, sex, region,
                                        symptoms, infected, vaccinated, vacc_inf),
              by = "unique_sample_name") %>%
    left_join(seropos$per_cov_prot, by = "unique_sample_name") %>%
    mutate(infected = as.character(infected),
           vaccinated = as.character(vaccinated)) %>%
    left_join(ifn_any_pos, by = "unique_sample_name") %>%
    select("Age group" = age_group, "Vaccine doses" = vaccinated, "Serocluster" = serocluster,
           sex, region, symptoms, infected, `N+`, `S+`, ifn_any_pos) %>%
    # Rename columns for table
    rename_with(\(x) {str_replace_all(x, "_", " ") %>% str_to_sentence()})
  
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
  
  tbl1 <- table1(~ `N+` + `S+` + `Age group` + Sex + Region +
                   Symptoms + `Vaccine doses` + Infected + `Ifn any pos` | Serocluster,
                 data = d, overall = F, extra.col = list(`P-value` = pvalue))
  
  return(tbl1)
}

# Heatmap of clusters
make_seroclust_heatmap <- function(seroclust, clust_cols) {
  mfi <- seroclust$data %>% rename_with(rename_serol)
  sinfo <- seroclust$serocluster
  
  # Viridis palette for heatmap cells
  col_fun <- colorRamp2(breaks = c(min(mfi[, -1]), max(mfi[, -1])), hcl_palette = "Viridis")
  
  # Sample annotation
  col_annot <- sinfo %$%
    columnAnnotation("Questionnaire" = vacc_inf,
                     # Define colours for each annotation
                     col = list("Questionnaire" = clust_cols %>% setNames(levels(vacc_inf))),
                     # Set annotation legends to be horizontal instead of vertical, increase font size
                     annotation_legend_param = list(list("nrow" = 2, labels_gp = gpar(fontsize = 16),
                                                         title_gp = gpar(fontsize = 18),
                                                         direction = "horizontal")) %>%
                       setNames(c("Questionnaire")),
                     annotation_name_gp = gpar(fontsize = 0),
                     annotation_name_side = "left")
  # Column split
  col_spl_hcl <- sinfo[match(mfi$unique_sample_name, sinfo$unique_sample_name), ] %>% pull(serocluster)
  col_spl_qst <- sinfo[match(mfi$unique_sample_name, sinfo$unique_sample_name), ] %>% pull(vacc_inf)
  
  by_hcl <- Heatmap(mfi %>% select(-unique_sample_name) %>% t(),
                    col = col_fun, show_column_names = F, name = "Z-score",
                    cluster_rows = T, clustering_distance_rows = "manhattan",
                    row_names_side = "left", row_dend_side = "right",
                    clustering_method_rows = "ward.D2",
                    top_annotation = col_annot, column_split = col_spl_hcl, cluster_column_slices = F,
                    heatmap_legend_param = list(legend_direction = "horizontal",
                                                title_gp = gpar(fontsize = 18),
                                                labels_gp = gpar(fontsize = 16)),
                    column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 20))
  by_qst <- Heatmap(mfi %>% select(-unique_sample_name) %>% t(),
                    col = col_fun, show_column_names = F, name = "Z-score",
                    cluster_rows = T, clustering_distance_rows = "manhattan",
                    row_names_side = "left", row_dend_side = "right",
                    clustering_method_rows = "ward.D2",
                    top_annotation = col_annot, column_split = col_spl_qst, cluster_column_slices = F,
                    heatmap_legend_param = list(legend_direction = "horizontal",
                                                title_gp = gpar(fontsize = 18),
                                                labels_gp = gpar(fontsize = 16)),
                    column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 20))
  
  return(list("hcl" = by_hcl, "qst" = by_qst))
}

# Boxplots of Abs in clusters, compared to questionnaire variables
make_seroclust_boxplt <- function(seroclust) {
  mfi <- seroclust$data
  sinfo <- seroclust$serocluster
  
  boxplt_dat <- mfi %>%
    left_join(sinfo %>% select(unique_sample_name, vacc_inf, serocluster),
              by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, vacc_inf, serocluster),
                 names_to = "ab", values_to = "zscore") %>%
    pivot_longer(cols = -c(unique_sample_name, ab, zscore), names_to = "name", values_to = "cluster") %>%
    filter(!is.na(cluster)) %>%
    # Order clusters as factor levels were lost when pivoting
    mutate(cluster = factor(cluster, levels = c(1:4, levels(sinfo$vacc_inf))))
  
  # Horizontal line at median of cluster 1 (V- I-)
  boxplt_hline <- boxplt_dat %>%
    filter(cluster %in% c("1", "V- I-")) %>%
    summarise(md = median(zscore), .by = c(name, ab))
  
  # Kruskal-Wallis test p-values and positions
  kwp <- boxplt_dat %>%
    nest(.by = c(ab, name)) %>%
    mutate(p = map(data, \(x) {
      kruskal.test(zscore ~ cluster, data = x)$p.value %>%
        signif(3)
    })) %>%
    mutate(ypos = map(data, \(x) {
      max(x$zscore) + 0.1 * diff(range(x$zscore))
    })) %>%
    select(-data) %>%
    unnest(c(p, ypos)) %>%
    mutate(xpos = 2.5)
  
  clust_boxplt <- ggplot(boxplt_dat, aes(x = cluster, y = zscore)) +
    geom_hline(data = boxplt_hline, mapping = aes(yintercept = md),
               colour = "grey", linetype = "dashed") +
    geom_boxplot() +
    geom_text(data = kwp, mapping = aes(x = xpos, y = ypos, label = p)) +
    # geom_boxplot(aes(fill = cluster), show.legend = F) +
    facet_grid2(rows = vars(ab), cols = vars(name), scales = "free", drop = T, independent = T) +
    # scale_fill_manual(values = c(clust_cols, clust_cols) %>%
    #                     setNames(c(1:4, levels(sinfo$vacc_inf)))) +
    theme_classic()
}

# Response and time for Abs in clusters
make_seroclust_responseplt <- function(seroclust, serol_dat, ab_cols) {
  
  # Use scaled MFI (Z-scores) as we are combining binders targeting the same protein for the trend lines
  time_plt_dat <- seroclust$data %>%
    rename_with(rename_serol) %>%
    as.data.frame() %>%
    # Cluster data
    left_join(seroclust$serocluster %>%
                select(unique_sample_name, serocluster),
              by = "unique_sample_name") %>%
  # Questionnaire data
    left_join(serol_dat$sinfo %>% select(unique_sample_name, vacc_inf, vacc_diff, inf_diff),
              by = "unique_sample_name") %>%
    # Add label for vaccination/infection status of majority of clusters (represented by anti-N and -S)
    mutate(clust_ns = case_when(serocluster == 1 ~ "N-S-",
                                serocluster == 2 ~ "N-S+",
                                serocluster == 3 ~ "N+S-",
                                serocluster == 4 ~ "N+S+")) %>%
    mutate(serocluster = paste0("Serocluster ", serocluster, " (", clust_ns, ")")) %>%
    select(-clust_ns) %>%
    # Reverse plotting order for figure
    mutate(serocluster = factor(serocluster, levels = rev(unique(as.character(serocluster))))) %>%
    # Make longer for Abs
    pivot_longer(cols = -c(unique_sample_name, vacc_inf, serocluster, vacc_diff, inf_diff),
                 names_to = "ab", values_to = "zscore") %>%
    # Make longer again with respect to variables to group by (to nest by for x axis)
    pivot_longer(cols = c(inf_diff, vacc_diff)) %>%
    filter(!is.na(value)) %>%
    # Add grouping variable for Abs (merge anti-N with each other, anti-S with each other for loess curves)
    mutate(ab_grp = case_when(str_detect(ab, "-N") ~ "Anti-N",
                              str_detect(ab, "-S|-R") ~ "Anti-S")) 
  
  # Antibody levels over time in selected clusters
  time_plt_hcl <- time_plt_dat %>%
    group_by(name) %>% nest() %>%
    mutate(plt = map(data, \(plt_dat) {
      plt_dat %>%
        filter(name == "inf_diff" & str_detect(serocluster, "4") |
               name == "vacc_diff" & str_detect(serocluster, "2")) %>%
        ggplot(aes(x = value, y = zscore)) +
        geom_point(aes(colour = ab, size = ab)) +
        scale_colour_manual(values = ab_cols) +
        scale_size_manual(values = c(
          "Anti-Na" = 2.5, "Anti-RBD" = 2.5, "Anti-Nc" = 3.5, "Anti-S1" = 3.5, "Anti-S1S2" = 3.5
        )) +
        labs(x = "Antibody", colour = "Antibody", size = "Antibody") +
        new_scale_colour() +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        geom_smooth(aes(colour = ab_grp), formula = y ~ x, method = "loess", se = T,
                    span = 1, alpha = 0.1, linewidth = 1.75) +
        scale_colour_manual(values = ab_cols) +
        facet_wrap2(~ serocluster, drop = T, axes = "all") +
        labs(x = switch(name, "inf_diff" = "Months since infection",
                        "vacc_diff" = "Months since vaccination"),
             y = "Z-score", colour = "Trend") +
        theme_classic(20)
    })) %>% pull(plt) %>% wrap_plots(guides = "collect")
  
  return(list("hierarchical" = time_plt_hcl))
}

# Look at samples that are infected or vaccinated but cluster in seroclusters that don't match their self-reported status
get_seroclust_mismatch <- function(seroclust, serol_in) {
  # Use non-scaled MFI
  mfi <- serol_in$mfi %>% filter(unique_sample_name %in% seroclust$data$unique_sample_name) %>%
    get_cov()
  
  # Table of mismatches
  tbl <- seroclust$serocluster %>%
    filter(vacc_inf != serocluster_vi) %>%
    count(serocluster, serocluster_vi, vacc_inf) %>%
    add_row(serocluster = "Total", n = sum(.$n))
  
  # Mark "mismatched" samples
  sinfo <- seroclust$serocluster %>%
    left_join(serol_in$sinfo %>% select(unique_sample_name, vaccine_company, symptoms,
                                        infected, vaccinated, vacc_diff, inf_diff),
              by = "unique_sample_name") %>%
    mutate(inf_mismatch = case_when(infected > 0 & str_detect(serocluster_vi, "I-") ~ "Mismatch", T ~ "Match"),
           vacc_mismatch = case_when(vaccinated > 0 & str_detect(serocluster_vi, "V-") ~ "Mismatch", T ~ "Match"))
  
  # Boxplots of anti-SARS-CoV-2 Abs for mismatched and not mismatched samples
  boxplot_dat <- mfi %>%
    left_join(sinfo %>% select(unique_sample_name, contains("mismatch"), vaccinated, infected),
              by = "unique_sample_name") %>%
    pivot_longer(cols = -c(unique_sample_name, contains("mismatch"), vaccinated, infected)) %>%
    mutate(name = rename_serol(name))
  
  boxplot1 <- boxplot_dat %>%
    select(inf_mismatch, infected, name, value) %>%
    filter(infected > 0) %>%
    ggplot(aes(x = inf_mismatch, y = value)) +
    geom_beeswarm(alpha = 0.5, colour = "grey50") +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    facet_wrap(~ name, nrow = 1) +
    labs(x = "Serocluster match, infected individuals", y = "Ab level [AU]") +
    theme_classic()
  
  boxplot2 <- boxplot_dat %>%
    select(vacc_mismatch, vaccinated, name, value) %>%
    filter(vaccinated > 0) %>%
    ggplot(aes(x = vacc_mismatch, y = value)) +
    geom_beeswarm(alpha = 0.5, colour = "grey50") +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    facet_wrap(~ name, nrow = 1) +
    labs(x = "Serocluster match, vaccinated individuals", y = "Ab level [AU]") +
    theme_classic()
  
  boxplot_combined <- boxplot1 / boxplot2
  
  # Barplots of number of doses and time between vaccination and sampling (average number of months)
  barplots <- sinfo %>%
    filter(vaccinated > 0) %>%
    select(vacc_mismatch, vaccinated, vacc_diff) %>%
    # All variables are categorical
    mutate(across(where(is.numeric), as.character)) %>%
    pivot_longer(cols = -vacc_mismatch) %>%
    # Order plots
    mutate(name = factor(name, levels = c("vaccinated", "vacc_diff"))) %>%
    group_by(name) %>% nest() %>% arrange(name) %>%
    mutate(plt = map(data, \(x) {
      comp <- comp_groups(x$vacc_mismatch, x[, "value", drop = T], discr_palette = "BuGn")$plt +
        labs(x = NULL,
             fill = switch(
               as.character(name), "vaccinated" = "Vaccine doses",
               "vacc_diff" = "Months since\nvaccination"
             )
        ) +
        coord_flip(clip = "off") +
        theme_classic(20)
      
      # Remove axis title and labels in the top plot
      if (name == "vaccinated") {
        comp <- comp + theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())
      }
      
      return(comp)
    })) %>%
    pull(plt) %>%
    wrap_plots(ncol = 1)
  
  return(list("boxplots" = boxplot_combined, "barplots" = barplots))
}

# Make a figure summarising serology analyses
make_serology_figure <- function(serol_boxplots, clust_heatmap, mismatch, serol_over_time) {
  
  # Boxplots of selected Abs
  abs_to_plot <- c(
    "Anti-S1S2 [AU]" = "anti_s1s2",
    "Anti-Nc [AU]" = "anti_nc",
    "Anti-RBD [AU]" = "anti_rbd",
    "Anti-EBNA1 [AU]" = "EBNA1"
  )
  
  selected_boxplots <- 
    # Collect four of the boxplots
    (imap(abs_to_plot, \(value, name) {
      serol_boxplots[[value]] + labs(x = NULL, y = name, fill = "Region") +
        theme_classic(20)
    }) %>%
      wrap_plots(., nrow = 2, guides = "collect"))
  
  # Grab the heatmap to be able to combine with ggplot2 objects
  clust_heatmap <- grid::grid.grabExpr(draw(clust_heatmap$hcl, merge_legend = T,
                                            heatmap_legend_side = "right",
                                            annotation_legend_side = "right"))
  
  # Combined plot
  plot_combined <- wrap_plots(
    selected_boxplots,
    clust_heatmap,
    mismatch$barplots & guides(fill = guide_legend(nrow = 2)),
    serol_over_time$hierarchical,
    nrow = 4, heights = c(1.2, 1, 0.8, 0.8)
  )
  
  # Save combined plot and return path as output
  plt_path <- make_fn("multianalyte_serology", "serology_figure.pdf")
  ggsave(plt_path, plot_combined, height = 18, width = 12)
  
  return(plt_path)
}

# Classification/prediction using serology
serol_classification <- function(serol_dat, seroclust, use_ifn = T, clust_cols) {
  d <- serol_dat$mfi %>%
    left_join(serol_dat$sinfo, by = "unique_sample_name") %>%
    filter(class == "Sample") %>%
    select(unique_sample_name, matches("anti_(s|rbd|n)"), contains("IFN")) %>%
    left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster), by = "unique_sample_name") %>%
    left_join(serol_dat$sinfo %>% select(unique_sample_name, age_group, sex, region, infected, vaccinated, vacc_inf), by = "unique_sample_name") %>%
    mutate(infected = as.character(infected), vaccinated = case_when(vaccinated > 0 ~ "1", T ~ "0"))
  
  resp <- c("serocluster", "age_group", "sex", "region", "infected", "vaccinated", "vacc_inf")
  pred <- paste0("anti_(s|rbd|n", ifelse(use_ifn, "|IFN", ""), ")")
  
  lasso_res <- sapply(resp, \(y) {
    # Omit NA
    d_in <- d %>% select(matches(pred), !!y) %>% filter(complete.cases(.))
    
    tidy_classification(x = d_in, response = !!y)
  }, simplify = F, USE.NAMES = T) %>%
    # Keep only the necessary results objects
    lapply(\(x) {x[c("metrics", "varimp", "roc_curve", "vip", "conf_plot")]})
  
  lasso_res$serocluster$roc_curve <- lasso_res$serocluster$roc_curve +
    scale_colour_manual(values = clust_cols)
  
  return(lasso_res)
}

