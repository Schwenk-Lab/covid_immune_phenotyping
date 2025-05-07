# Data wrangling target functions

# Read raw serology data
read_serol_raw <- function(serol_files, s_low_conc) {
  d <- list(
    mfi = read_csv(str_subset(serol_files, "mfi"), show_col_types = F),
    sinfo = read_csv(str_subset(serol_files, "sinfo"), show_col_types = F),
    binfo = read_csv(str_subset(serol_files, "binfo"), show_col_types = F)
  )
  
  # Fix sample names if any sample name is missing a leading zero, causing mismatch with other data types
  short_name <- d$sinfo %>% filter(class == "Sample") %>%
    mutate(name_len = nchar(sample_name)) %>%
    filter(name_len == 2) %>%
    pull(sample_name)

  # Add a leading zero
  d$mfi <- d$mfi %>%
    mutate(unique_sample_name = case_when(unique_sample_name %in% paste0("X", short_name) ~ paste0("X0", str_remove(unique_sample_name, "X")), T ~ unique_sample_name))
  rownames(d$mfi) <- d$mfi$unique_sample_name

  d$sinfo <- d$sinfo %>% mutate(unique_sample_name = case_when(sample_name %in% short_name ~ paste0("X0", sample_name), T ~ unique_sample_name),
                                sample_name = case_when(sample_name %in% short_name ~ paste0("0", sample_name), T ~ sample_name))
  rownames(d$sinfo) <- d$sinfo$unique_sample_name

  
  # Filter out analytes that will not be used and rename to be more readable
  # Pick out binders for serology
  if (s_low_conc) {
    b_incl <- c("unique_sample_name",
                str_subset(colnames(d$mfi),"(Spike|RBD).*0.4$"),
                str_subset(colnames(d$mfi), "Nucleo.*[[:alpha:]]$"),
                str_subset(colnames(d$mfi), "IFN"),
                str_subset(colnames(d$mfi), "EBNA|IgG|IgM|IgA"),
                str_subset(colnames(d$mfi), "Bare"))
  } else {
    b_incl <- c("unique_sample_name",
                str_subset(colnames(d$mfi),"(Spike).*[^4]$"),
                "RBD",
                str_subset(colnames(d$mfi),"(Nucleo).*[^4]$"),
                str_subset(colnames(d$mfi), "IFN"),
                str_subset(colnames(d$mfi), "EBNA|IgG|IgM|IgA"),
                str_subset(colnames(d$mfi), "Bare"))
  }
  
  d$binfo <- d$binfo %>%
    filter(unique_antigen_name %in% b_incl) %>%
    mutate(unique_antigen_name = case_when(str_detect(unique_antigen_name, "Nucleocapsid\\.C") ~ "anti_nc",
                                           str_detect(unique_antigen_name, "Nucleocapsid_Acro") ~ "anti_na",
                                           str_detect(unique_antigen_name, "Spike\\.S1S2\\.foldon") ~ "anti_s1s2",
                                           str_detect(unique_antigen_name, "Spike\\.S1") ~ "anti_s1",
                                           str_detect(unique_antigen_name, "RBD") ~ "anti_rbd",
                                           str_detect(unique_antigen_name, "IgM") ~ "anti_human_igm",
                                           str_detect(unique_antigen_name, "IgA") ~ "anti_human_iga",
                                           str_detect(unique_antigen_name, "human\\.IgG") ~ "anti_human_igg",
                                           str_detect(unique_antigen_name, "Rabbit\\.IgG") ~ "anti_rabbit_igm",
                                           T ~ unique_antigen_name))
  d$mfi <- d$mfi %>% select(!!b_incl) %>%
    rename_with(\(x) {
      case_when(str_detect(x, "Nucleocapsid\\.C") ~ "anti_nc",
                str_detect(x, "Nucleocapsid_Acro") ~ "anti_na",
                str_detect(x, "Spike\\.S1S2\\.foldon") ~ "anti_s1s2",
                str_detect(x, "Spike\\.S1") ~ "anti_s1",
                str_detect(x, "RBD") ~ "anti_rbd",
                str_detect(x, "IgM") ~ "anti_human_igm",
                str_detect(x, "IgA") ~ "anti_human_iga",
                str_detect(x, "human\\.IgG") ~ "anti_human_igg",
                str_detect(x, "Rabbit\\.IgG") ~ "anti_rabbit_igm",
                T ~ x)
    })
  
  d$comment <- "Raw MFI data for both 384-plates, DBS vaccination serology"
  
  return(d[c("mfi", "sinfo", "binfo", "comment")])
}

# Normalise serology data
norm_serol <- function(dat_in, script) {
  
  # Input data for normalisation script
  write_csv(dat_in$mfi %>%
              mutate(Sample = case_when(grepl("^X", unique_sample_name) ~ "Sample",
                                        grepl("NL\\.neg", unique_sample_name) ~ "Negctrl",
                                        T ~ "Empty")) %>%
              relocate(Sample),
            file = "../data/mixmod_in.csv")
  
  # Run normalisation in environment with julia installed
  system(paste0("conda run -n julia_env julia ", script))
  
  # Load normalised data, filter out old columns
  norm_dat <- read_csv("../data/WithAdjusted.csv", show_col_types = F) %>%
    select(unique_sample_name, ends_with("_corrected_log2")) %>%
    rename_with(\(x) {
      str_replace(x, "_corrected_log2", "")
    })
  
  dat_out <- dat_in
  dat_out$mfi <- norm_dat
  dat_out$comment <- "Serology MFI normalised for bare-bead correlation with mixed model"
  
  return(dat_out)
}

# Plot raw vs normalised serology
plot_serol_norm <- function(dat_in) {
  # Take dat as input to make the step dependent on the normalisation step, but the data is loaded from the file to get both old and new values
  
  norm_dat <- read_csv("../data/WithAdjusted.csv", show_col_types = F)
  lines <- read_csv("../data/line_df.csv", show_col_types = F)
  
  nonblank_ind <- norm_dat$Sample != "Empty"
  prots <- str_subset(colnames(norm_dat)[-c(1, 2)], "_corrected_log2", negate = T)
  
  plt_file <- make_fn("serology_normalisation", "norm_before_after.pdf")
  pdf(plt_file, onefile = T, width = 12, height = 5)
  
  for (i in prots) {
    
    plt_dat <- data.frame("samp_name" = norm_dat$unique_sample_name) %>%
      mutate("Raw" = log2(norm_dat[, i, drop = T]),
             "Normalised" = norm_dat[, paste0(i, "_corrected_log2"), drop = T],
             "Bare.bead.2" = log2(norm_dat[, "Bare.bead.2", drop = T])) %>%
      mutate(Sample_type = case_when(grepl("blank", samp_name, ignore.case = T) ~ "Blank",
                                     grepl("pool", samp_name, ignore.case = T) ~ "Pool",
                                     grepl("NL\\.neg", samp_name, ignore.case = T) ~ "Neg ctrl",
                                     grepl("NL\\.pos", samp_name, ignore.case = T) ~ "Pos ctrl",
                                     T ~ "Sample"))
    
    p_before <- ggplot(plt_dat) +
      geom_point(aes(x = Bare.bead.2, y = Raw, colour = Sample_type)) +
      geom_segment(data = lines %>% filter(protein == i),
                   aes(x = mi_val, xend = ma_val,
                       y = o_m_val * mi_val + o_c_val, yend = o_m_val * ma_val + o_c_val)) +
      scale_colour_brewer(palette = "Set2")
    
    # Blank deviation for plot subtitle
    mean_blank_dev <- lines %>% filter(protein == i) %>%
      summarise(mbd = mean(plt_dat$Raw[!nonblank_ind] - (o_m_val * plt_dat$Bare.bead.2[!nonblank_ind] + o_c_val))) %>%
      pull(mbd)
    
    p_after <- ggplot(plt_dat) +
      geom_point(aes(x = Bare.bead.2, y = Normalised, colour = Sample_type)) +
      scale_colour_brewer(palette = "Set2")
    
    p_both <- p_before + p_after +
      plot_layout(guides = "collect") +
      plot_annotation(title = i, subtitle = paste0(
        "m=", lines[lines$protein == i, "o_m_val"], "; Mean Blank deviation=", mean_blank_dev
      )) &
      theme_bw() &
      theme(axis.line = element_line(linewidth = 1),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 15))
    
    print(p_both)
  }
  dev.off()
  
  return(plt_file)
}

# Read Olink data and perform bridge normalisation
bridge_normalisation_olink <- function(olink_file1, olink_file2, olink_binfo) {
  
  # Data from two different Olink runs
  d1 <- list(sinfo = read_csv(str_subset(olink_file1, "sinfo"), show_col_types = F),
             long_data = read_csv(str_subset(olink_file1, "long_data"), show_col_types = F))
  d2 <- list(sinfo = read_csv(str_subset(olink_file2, "sinfo"), show_col_types = F),
             long_data = read_csv(str_subset(olink_file2, "long_data"), show_col_types = F))
  
  # Perform bridge normalisation
  # Get overlapping binders, skip failed binders with no values (missing frequency is NA)
  bridge_binders <- intersect(
    d1$long_data %>% filter(!is.na(missing_freq), panel == "Inflammation") %>%
      distinct(assay) %>% pull(assay),
    d2$long_data %>% filter(!is.na(missingfreq)) %>% distinct(assay) %>% pull(assay)
  )
  
  bridge_samples <- intersect(d1$sinfo$unique_sample_name, d2$sinfo$unique_sample_name)
  
  bridge_diff <- rbind(
    d1$long_data %>% filter(unique_sample_name %in% bridge_samples,
                            panel == "Inflammation", assay %in% bridge_binders) %>%
      dplyr::select(unique_sample_name, assay, npx) %>% mutate(proj = "1"),
    d2$long_data %>% filter(matchid %in% bridge_samples, assay %in% bridge_binders) %>%
      dplyr::select(matchid, assay, npx) %>%
      rename(unique_sample_name = matchid) %>% mutate(proj = "2")
  ) %>%
    # Per sample difference for each protein
    group_by(unique_sample_name, assay) %>%
    summarise(prot_diff = npx[proj == "2"] - npx[proj == "1"], .groups = "keep") %>%
    # Per protein median of differences (normalisation factor)
    ungroup(unique_sample_name) %>%
    summarise(med_diff = median(prot_diff))
  
  # Add the normalisation factor to the Olink subset dataset
  bn_long <- d1$long_data %>% filter(panel == "Inflammation",
                                     assay %in% bridge_binders) %>%
    left_join(d1$sinfo %>% dplyr::select(unique_sample_name, sample_type), by = "unique_sample_name")
  bn_long$npx <- bn_long$npx + bridge_diff[match(bn_long$assay, bridge_diff$assay), "med_diff", drop = T]
  # Also update LOD
  bn_long$lod <- bn_long$lod + bridge_diff[match(bn_long$assay, bridge_diff$assay), "med_diff", drop = T]
  
  # Bind together the two datasets
  col_select <- c("unique_sample_name", "index", "olink_id", "uniprot", "assay", "missing_freq", "panel", "panel_lot_nr", "plate_id", "qc_warning", "lod", "npx", "normalization", "assay_warning")
  # Pooled samples
  bn_pools <- rbind(
    d2$long_data %>% filter(sample_type == "POOL", assay %in% bridge_binders) %>%
      rename(unique_sample_name = matchid, olink_id = olinkid, missing_freq = missingfreq, plate_id = plateid) %>%
      dplyr::select(!!col_select),
    bn_long %>%
      filter(!unique_sample_name %in% bridge_samples, str_detect(sample_type, "pool")) %>%
      dplyr::select(!!col_select)
  )
  # Normal samples
  bn_combine <- rbind(
    d2$long_data %>% filter(sample_type == "SAMPLE", assay %in% bridge_binders) %>%
      rename(unique_sample_name = matchid, olink_id = olinkid, missing_freq = missingfreq, plate_id = plateid) %>%
      dplyr::select(!!col_select),
    bn_long %>%
      filter(!unique_sample_name %in% bridge_samples, sample_type == "sample") %>%
      dplyr::select(!!col_select)
  )
  
  return(list(
    "bridge_samples" = bridge_samples,
    "normalisation_factors" = bridge_diff,
    "samples" = bn_combine,
    "pools" = bn_pools,
    "binfo" = read_csv(olink_binfo, show_col_types = F)
  ))
}

# Make data list from bridge normalisation results
get_olink_raw <- function(olink_bn, serol_in) {
  bn_npx <- olink_bn$samples %>% select(unique_sample_name, assay, npx) %>%
    pivot_wider(id_cols = unique_sample_name, names_from = assay, values_from = npx)
  bn_npx_pools <- olink_bn$pools %>% select(unique_sample_name, assay, npx) %>%
    pivot_wider(id_cols = unique_sample_name, names_from = assay, values_from = npx)
  # Recalculate missing frequencies
  mf <- olink_bn$samples %>%
    summarise(missingfreq = sum(npx < lod) / n(), .by = assay)
  
  sinfo_out <- serol_in$sinfo %>% filter(unique_sample_name %in% bn_npx$unique_sample_name)
  binfo_out <- olink_bn$binfo %>% select(-lod, -missingfreq) %>%
    left_join(mf, by = "assay") %>%
    filter(!is.na(missingfreq))
  
  dat <- list(
    "npm" = bn_npx,
    "npm_pool" = bn_npx_pools,
    "sinfo" = sinfo_out,
    "binfo" = binfo_out %>% filter(assay %in% colnames(bn_npx)),
    "comment" = "Bridge normalised Olink data."
  )
  
  return(dat)
}

# Filter Olink data
filter_olink <- function(prot_in) {
  prot_in$binfo <- prot_in$binfo %>% filter(missingfreq < 0.5)
  prot_in$npm <- prot_in$npm %>% select(unique_sample_name, !!prot_in$binfo$assay)
  # Exclude pooled samples
  prot_in <- prot_in[-which(names(prot_in) == "npm_pool")]
  
  return(prot_in)
}

# Read Alamar data
read_alamar_raw <- function(alamar_files) {
  list(
    npm = read_csv(str_subset(alamar_files, "npq"), show_col_types = F),
    sinfo = read_csv(str_subset(alamar_files, "sinfo"), show_col_types = F),
    binfo = read_csv(str_subset(alamar_files, "binfo"), show_col_types = F),
    lod = read_csv(str_subset(alamar_files, "lod"), show_col_types = F),
    comment = "Alamar NULISA raw data"
  )
}

# Filter Alamar data
filter_alamar <- function(prot_in) {
  # Keep proteins with >50% of binders above LOD
  prot_in$binfo <- prot_in$binfo %>% filter(overall > 50)
  # Exclude pools, controls
  prot_in$npm <- prot_in$npm %>%
    select(unique_sample_name, !!prot_in$binfo$assay) %>%
    filter(str_detect(unique_sample_name, "^X"))
  prot_in$lod <- prot_in$lod %>% filter(assay %in% prot_in$binfo$assay)
  prot_in$sinfo <- prot_in$sinfo %>%
    filter(unique_sample_name %in% prot_in$npm$unique_sample_name)
  
  return(prot_in)
}

# Apply ProtPQN
apply_norm_protpqn <- function(dat_in) {
  norm_dat <- dat_in$npm %>% column_to_rownames("unique_sample_name") %>%
    2 ** . %>% as.matrix() %>% protpqn() %>% as.data.frame() %>% log2() %>%
    rownames_to_column("unique_sample_name")
  
  dat_out <- dat_in
  dat_out$npm <- norm_dat
  
  return(dat_out)
}

# Outlier detection by sample median and IQR (similar to OlinkAnalyze::olink_qc_plot())
med_iqr_plt <- function(dat_in, pt_size = 2, threshold_sds = 3, lab_outl = T) {
  # dat_in, long format data with sample names ('unique_sample_name'), protein measurements ('npm'), and protein panel names ('panel')
  # Protein names are not necessary but could help for checking that the input is correct
  # pt_size, numeric for point size
  # lab_outl, boolean, T to display sample names of outliers
  
  # Compute medians and IQRs of samples, per protein panel
  iqr_median_dat <- dat_in %>%
    group_by(unique_sample_name, panel) %>%
    summarise(iqr = IQR(npm), median = median(npm), .groups = "keep") %>%
    ungroup()
  
  # Thresholds for outliers
  iqr_median_thr <- iqr_median_dat %>%
    group_by(panel) %>%
    summarise(upper_iqr = mean(iqr) + threshold_sds * sd(iqr),
              lower_iqr = mean(iqr) - threshold_sds * sd(iqr),
              upper_median = mean(median) + threshold_sds * sd(median),
              lower_median = mean(median) - threshold_sds * sd(median),
              .groups = "keep")
  
  # Add column saying whether a sample is an outlier or not
  iqr_median_dat <- left_join(iqr_median_dat, iqr_median_thr, by = c("panel")) %>%
    mutate(outlier = case_when(iqr > upper_iqr | iqr < lower_iqr | median > upper_median | median < lower_median ~ "outlier"))
  
  # Turn thresholds table into long format for plotting
  iqr_median_thr <- pivot_longer(iqr_median_thr, -c(panel), names_to = "threshold", values_to = "value")
  
  # Make plot, nest by protein panel to keep same style as the other plots
  med_iqr_plt <- iqr_median_dat %>% group_by(panel) %>% nest() %>%
    mutate(plt = map(data, ~ {
      # Get panel name to avoid conflicts arising from using the same name on variables
      cur_panel <- panel
      
      # Make data frame for samples to draw lines between
      replicate_ids <- .x %>%
        filter(str_detect(unique_sample_name, "_[:alpha:]$")) %>%
        mutate(unique_sample_name = str_remove(unique_sample_name, "_[:alpha:]$")) %>%
        pull(unique_sample_name)
      
      # Deal with case of no replicates
      if (length(replicate_ids) > 0) {
        line_df <- .x %>%
          filter(str_detect(unique_sample_name, paste0(replicate_ids, collapse = "|"))) %>%
          mutate(line_id = str_remove(unique_sample_name, "_[:alpha:]$"))
      } else {
        line_df <- .x %>% mutate(line_id = "") %>% slice(0)
      }
      
      thr_df <- iqr_median_thr %>% filter(panel == cur_panel)
      
      p <- .x %>% mutate(panel = cur_panel) %>%
        ggplot(aes(x = median, y = iqr)) +
        geom_line(aes(group = line_id), colour = "grey", data = line_df) +
        geom_point(size = pt_size) +
        geom_hline(aes(yintercept = value), data = thr_df %>% filter(str_detect(threshold, "iqr")),
                   colour = "grey75", linetype = "dashed") +
        geom_vline(aes(xintercept = value), data = thr_df %>% filter(str_detect(threshold, "median")),
                   colour = "grey75", linetype = "dashed") +
        geom_smooth(data = .x %>% filter(is.na(outlier)), formula = y ~ x, method = "lm", alpha = 0.3) +
        labs(x = "Median", y = "IQR") +
        facet_wrap(~ panel) +
        theme_classic()
      
      # Label outliers
      if (lab_outl) {p <- p + geom_text_repel(aes(label = unique_sample_name), data = .x %>% filter(outlier == "outlier"), size = 3)}
      
      return(p)
    })) %>%
    pull(plt)
  
  return(list("plots" = med_iqr_plt, "data" = iqr_median_dat))
}

# Remove samples detected as outliers
outlier_exclude <- function(dat_in, outl_det) {
  dat_out <- dat_in
  outliers <- outl_det$data %>% filter(!is.na(outlier)) %>%
    pull(unique_sample_name)
  dat_out$npm <- dat_out$npm %>% filter(!unique_sample_name %in% outliers)
  
  return(dat_out)
}

# Make full proteomics data set with all proteins
combine_proteomics <- function(olink_in, alamar_in) {
  
  # Combine replicates if enough volume
  alm_comb <- combine_repl(alamar_in)
  
  # Check overlapping proteins, give them an alternative name
  ol_prot <- intersect(colnames(olink_in$npm)[-1],
                       colnames(alm_comb$npm)[-1])
  binfo_o <- olink_in$binfo %>%
    mutate(unique_protein_name = case_when(assay %in% ol_prot ~ paste0(assay, "_1"), T ~ assay))
  binfo_a <- alm_comb$binfo %>%
    mutate(unique_protein_name = case_when(assay %in% ol_prot ~ paste0(assay, "_2"), T ~ assay))
  
  # Merge
  npm_both <- olink_in$npm %>%
    # Update column names
    rename_with(\(x) {
      c("unique_sample_name",
        binfo_o[match(x[-1], binfo_o$assay), ] %>% pull(unique_protein_name))
    }) %>%
    left_join(alm_comb$npm %>%
                rename_with(\(x) {
                  c("unique_sample_name",
                    binfo_a[match(x[-1], binfo_a$assay), ] %>% pull(unique_protein_name))
                }), by = "unique_sample_name") %>%
    filter(complete.cases(.))
  
  # Protein info list with all proteins, using only common columns
  binfo_both <- rbind(
    binfo_o %>% select(unique_protein_name, assay, uniprot, panel, assay_id = olinkid) %>% mutate(method = "olink"),
    binfo_a %>% select(unique_protein_name, assay, uniprot, panel, assay_id = alamartargetid) %>% mutate(method = "alamar")
  )
  
  dat_out <- list(npm = npm_both,
                  sinfo = olink_in$sinfo %>% filter(unique_sample_name %in% npm_both$unique_sample_name),
                  binfo_a = binfo_a, binfo_o = binfo_o, binfo = binfo_both,
                  comment = "Total proteomics data for 500 proteins from Olink and Alamar, with only high-detectability (>50%) proteins, normalised with ProtPQN, and with outliers removed.")
  
  return(dat_out)
}

# Adjust proteomics data for confounding variables
adjust_proteomics <- function(prot_in) {
  prot_out <- prot_in
  npm <- prot_out$npm

  # Adjust data for sex
  npm_out <- npm %>%
    left_join(prot_in$sinfo %>% select(unique_sample_name, sex), by = "unique_sample_name") %>%
    filter(!is.na(sex)) %>%
    pivot_longer(cols = -c(unique_sample_name, sex)) %>%
    nest(.by = name) %>%
    mutate(npm_adj = map(data, \(x) {
      lm(value ~ sex, data = x) %>% resid()
    })) %>%
    unnest(c(data, npm_adj)) %>%
    pivot_wider(id_cols = unique_sample_name, names_from = name, values_from = npm_adj)

  prot_out$npm <- npm_out
  prot_out$sinfo <- prot_in$sinfo %>% filter(unique_sample_name %in% npm_out$unique_sample_name)
  prot_out$comment <- "Total proteomics data for 500 proteins from Olink and Alamar, with only high-detectability (>50%) proteins, normalised with ProtPQN, with outliers removed, and adjusted for sex."

  return(prot_out)
}

# Tables of cohort demographics
make_demographics_tbl <- function(data_in) {
  # data_in, list containing data frame with sample info (e.g. serology data)
  
  d <- data_in$sinfo %>% filter(class == "Sample") %>%
    select(age_group, sex, region, infected, vaccinated, vacc_inf) %>%
    # Prepare variables for table
    mutate(infected = case_when(infected == 0 ~ "No", infected == 1 ~ "Yes"),
           vaccinated = as.character(vaccinated)) %>%
    # Rename columns for table
    rename(`Immune group` = vacc_inf, `Vaccine doses` = vaccinated) %>%
    rename_with(\(x) {str_replace(x, "_", " ") %>% str_to_sentence()})
  
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
  
  by_region <- table1(~ . | Region, data = d, overall = F, extra.col = list(`P-value` = pvalue))
  by_age <- table1(~ . | `Age group`, data = d, overall = F, extra.col = list(`P-value` = pvalue))
  
  return(list("by_region" = by_region, "by_age" = by_age))
}

# Barplot versions of cohort summary tables
make_demographics_bar <- function(data_in) {
  # data_in, list with sample info data frame (e.g. serology data)
  
  d <- data_in$sinfo %>% filter(class == "Sample")
  
  # Age and sex distributions in the regions
  age_sex_bar <- d %>% select(age_group, sex, region) %>%
    pivot_longer(cols = -region) %>%
    nest_by(name) %>%
    mutate(plt = list(
      data %>% filter(complete.cases(.)) %>%
        # Convert to percentages
        count(region, value) %>%
        mutate(percent = n / sum(n) * 100, .by = region) %>%
        ggplot(aes(x = value)) +
        geom_bar(aes(y = percent, fill = region), stat = "identity", position = "dodge") +
        scale_fill_brewer(palette = "Paired") +
        labs(x = str_to_sentence(str_replace(name, "_", " ")), y = "Percent [%]", fill = "Region") +
        theme_classic()
    )) %>%
    pull(plt)
  
  # Distributions of infected, vaccinated at least once, across age groups in each region
  inf_vacc_bar <- d %>% select(age_group, infected, vaccinated, region) %>%
    mutate(vaccinated = case_when(vaccinated > 0 ~ 1, T ~ vaccinated)) %>%
    summarise(infected = sum(infected, na.rm = T) / n() * 100,
              vaccinated = sum(vaccinated, na.rm = T) / n() * 100,
              .by = c(age_group, region)) %>%
    pivot_longer(cols = -c(age_group, region)) %>%
    nest_by(name) %>%
    mutate(plt = list(
      data %>% filter(!is.na(value)) %>%
        ggplot(aes(x = age_group, y = value, fill = region)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_brewer(palette = "Paired") +
        labs(x = "Age group", y = paste0(str_to_sentence(name), " [%]"), fill = "Region") +
        theme_classic()
    )) %>%
    pull(plt)
  
  return(list("age_sex" = age_sex_bar, "inf_vacc" = inf_vacc_bar))
  
}

# Barplots of sampling, infection, and vaccination timings
make_timings_bar <- function(data_in) {
  # data_in, list with sample info data frame (serology data)
  
  # Convert labels to be consistent between variables
  timings <- data_in$sinfo %>%
    filter(!is.na(consent_date)) %>%
    select(inf_date_latest, vacc_date_latest, consent_date) %>%
    mutate(inf_date_latest = case_when(!is.na(inf_date_latest) ~ str_replace_all(inf_date_latest, "-", "~") %>% str_replace_all(., "/", "-(") %>% paste0("20", ., ")"))) %>%
    mutate(vacc_date_latest = str_replace(vacc_date_latest, "/", "-")) %>%
    mutate(consent_date = as.Date(consent_date)) %>%
    mutate(inf_date_latest = case_when(inf_date_latest == "2021-(09)" ~ "2021-(07~09)", T ~ inf_date_latest))
  
  # Brewer Pastel2 colour scale for bar fills
  fill_cols <- RColorBrewer::brewer.pal(3, "Pastel2")
  
  # Sampling/consent date
  p1 <- ggplot(timings, aes(x = consent_date)) +
    geom_histogram(fill = fill_cols[1], binwidth = 30) +
    scale_x_date(date_breaks = "2 months", limits = c(as.Date("2020-01-01"), as.Date("2021-09-01")), date_labels = "%Y-%m") +
    labs(x = "Sampling date", y = "Count") +
    theme_classic()
  
  # Latest infection date
  p2 <- ggplot(timings %>% filter(!is.na(inf_date_latest)), aes(x = inf_date_latest)) +
    geom_histogram(stat = "count", fill = fill_cols[2]) +
    labs(x = "Date of latest infection", y = "Count") +
    theme_classic()
  
  # Latest vaccination date
  # To get the exact same x scales as the sampling plot, set the month intervals to be at the 15th of each month
  # then specify bins to contain 30 days, producing the same bars as would be produced in a barplot with
  # categorical x axis
  p3_dat <- timings %>% filter(!is.na(vacc_date_latest)) %>%
    mutate(vacc_date_latest = paste0(as.character(vacc_date_latest), "-15") %>% as.Date())
  p3 <- ggplot(p3_dat, aes(x = vacc_date_latest)) +
    geom_histogram(fill = fill_cols[3], binwidth = 30) +
    scale_x_date(date_breaks = "2 months",
                 limits = c(as.Date("2020-01-01"), as.Date("2021-09-01")), date_labels = "%Y-%m") +
    labs(x = "Latest vaccination date", y = "Count") +
    theme_classic()
  
  p1 / p2 / p3 & theme(panel.grid.major.y = element_line(colour = "grey90"))
}

# Look at distributions of anti-SARS-CoV-2 Ab positivity, anti-IFN AAb positivity, and infection, vaccination
# across age groups and sexes
check_age_sex_pos <- function(serol_dat, serol_seropos, ifn_seropos) {
  dat <- serol_dat$sinfo %>%
    # Get samples and their infection/vaccination statuses
    filter(class == "Sample") %>%
    select(unique_sample_name, infected, vaccinated) %>%
    mutate(infected = case_when(infected == 0 ~ "No", infected == 1 ~ "Yes"),
           vaccinated = case_when(vaccinated == 0 ~ "No", vaccinated %in% c(1, 2) ~ "Yes")) %>%
    # Add covid Ab seropositivity and recode
    left_join(serol_seropos$per_cov_prot, by = "unique_sample_name") %>%
    mutate(across(where(is.logical), \(x) {case_when(x ~ "Positive", T ~ "Negative")})) %>%
    # Convert to long format to add IFN AAb seropositivity (already in long format)
    pivot_longer(cols = -unique_sample_name, names_to = "var", values_to = "value") %>%
    bind_rows(ifn_seropos %>% select(unique_sample_name, var = name, value = seroclass) %>% filter(str_detect(var, "^IFN"))) %>%
    # Add age and sex information
    left_join(serol_dat$sinfo %>% select(unique_sample_name, age_group, sex), by = "unique_sample_name") %>%
    filter(!is.na(sex), !is.na(value))
  
  # Test differences in each age group
  ftest <- dat %>%
    nest(.by = c(var, age_group)) %>%
    mutate(test = map(data, \(x) {
      # NA if only one level in data
      if (length(unique(x$value)) == 1) return(list(p.value = NA, conf.int = NA,  estimate = NA))
      fisher.test(x$value, x$sex)
    })) %>%
    # Pick out p-value, confidence intervals and OR estimate and compute FDR-corrected p-values
    mutate(p_value = map(test, "p.value"),
           fdr = p.adjust(p_value, method = "fdr"),
           l_ci = map(test, \(x) x$conf.int[1]),
           u_ci = map(test, \(x) x$conf.int[2]),
           estimate = map(test, "estimate")) %>%
    select(-data, -test) %>%
    unnest(c(p_value, fdr, l_ci, u_ci, estimate))
  
  # Make plots, one with all IFN AAbs and N + S Abs, one with infection + vaccination
  # Split Ab plot into multiple due to large size
  plt_order <- dat %>% distinct(var) %>%
    filter(!var %in% c("infected", "vaccinated")) %>%
    # Order IFN + N + S alphabetically, numerically if same kind of IFN
    mutate(char = str_extract(var, "[:alpha:]+"),
           num = str_extract(var, "\\d+$") %>% as.numeric()) %>%
    arrange(char, num) %>%
    mutate(var = factor(var, levels = var)) %>%
    arrange(var) %>% select(var) %>%
    # Split plot into four due to size
    mutate(facet = ceiling(seq_along(var) / ceiling(length(var) / 4)))
  
  plt_abs <- dat %>%
    filter(!var %in% c("infected", "vaccinated")) %>%
    # Count for plotting
    count(var, value, age_group, sex) %>%
    mutate(percent = n / sum(n) * 100, .by = c(var, age_group, sex)) %>%
    # Order Abs
    mutate(var = factor(var, levels = levels(plt_order$var))) %>%
    arrange(var) %>%
    # Plotting order
    left_join(plt_order, by = "var") %>%
    group_by(facet) %>% nest() %>%
    mutate(plt = map(data, \(dat) {
      p <- ggplot(dat, aes(x = percent, y = age_group, fill = value)) +
        geom_col(position = "stack", colour = "white", width = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
        scale_fill_manual(values = c("Negative" = "grey", "Positive" = "#4876FF")) +
        # scale_fill_brewer(palette = "Set2", direction = -1) +
        # Show x axis on all inner plots but don't write out axis text
        facet_grid(rows = vars(var), cols = vars(sex), scales = "free_x",
                   axes = "all_x", axis.labels = "all_y", drop = T) +
        # Reverse scale of x axis for left facet, adjust axis limits so the facets both start at exactly 0
        facetted_pos_scales(x = list(scale_x_reverse(expand = expansion(mult = c(0.05, 0))), scale_x_continuous(expand = expansion(mult = c(0, 0.05))))) +
        labs(x = "Percent", y = "Age", fill = "Serostatus") +
        theme_classic(14) +
        # Reduce spacing between facets so the left and right bars touch
        theme(panel.spacing = unit(0, "pt"),
              strip.text = element_text(size = 14),
              strip.background = element_blank(),
              panel.grid.major.x = element_line())
      
      # Remove y axis label if inner plot
      if (facet != 1) {p <- p + theme(axis.text.y = element_blank())}
      
      return(p)
    })) %>% pull(plt) %>% wrap_plots(nrow = 1, guides = "collect") + plot_layout(axis_titles = "collect")
  
  plt_vaccinf <- dat %>%
    filter(var %in% c("infected", "vaccinated")) %>%
    # Count for plotting
    count(var, value, age_group, sex) %>%
    mutate(percent = n / sum(n) * 100, .by = c(var, age_group, sex)) %>%
    mutate(var = str_to_sentence(var)) %>%
    ggplot(aes(x = percent, y = age_group, fill = value)) +
    geom_col(position = "stack", colour = "white", width = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    scale_fill_manual(values = c("No" = "grey", "Yes" = "#4876FF")) +
    # scale_fill_brewer(palette = "Set2", direction = -1) +
    facet_grid(rows = vars(var), cols = vars(sex), scales = "free", axes = "all_x", axis.labels = "all_y") +
    facetted_pos_scales(x = list(scale_x_reverse(expand = expansion(mult = c(0.05, 0))), scale_x_continuous(expand = expansion(mult = c(0, 0.05))))) +
    labs(x = "Percent", y = "Age", fill = "Status") +
    theme_classic(14) +
    theme(panel.spacing = unit(0, "pt"),
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          panel.grid.major.x = element_line())
  
  return(list("test" = ftest, "plt_abs" = plt_abs, "plt_vaccinf" = plt_vaccinf))
}

# Make plots showing the data distribution of each measured analyte
plot_distrib_all <- function(serol_dat, prot_dat, serol_seropos, ifn_seropos) {
  
  # Combine serology and proteomics measurements into one table for plotting
  
  plt_dat_serol <- serol_dat$mfi %>%
    select(unique_sample_name, contains(c("anti_n", "anti_s", "rbd", "IFN"))) %>%
    filter(unique_sample_name %in% (serol_dat$sinfo %>%
                                      filter(class == "Sample") %>%
                                      pull(unique_sample_name))) %>%
    pivot_longer(cols = -unique_sample_name) %>%
    # Add info on serostatus
    left_join(serol_seropos$per_ab %>% select(name, unique_sample_name, seroclass_serol = seroclass),
              by = c("name", "unique_sample_name")) %>%
    left_join(ifn_seropos %>% select(name, unique_sample_name, seroclass_ifn = seroclass),
              by = c("name", "unique_sample_name")) %>%
    unite("seroclass", seroclass_serol, seroclass_ifn, sep = "", na.rm = T) %>%
    mutate(type = case_when(str_detect(name, "IFN") ~ "ifn", T ~ "cov")) %>%
    mutate(name = rename_serol(name)) %>%
    arrange(type, name)
  
  plt_dat_prot <- prot_dat$npm %>%
    pivot_longer(cols = -unique_sample_name) %>%
    mutate(seroclass = NA, type = "prot") %>%
    arrange(name)
  
  # Make a plot per binder, colour by seropositivity
  plt <- rbind(plt_dat_serol, plt_dat_prot) %>%
    # Some interferons are included in both AAb assays and proteomics assays, scale by both name and assay type
    # mutate(value = scale(value), .by = c(name, type)) %>%
    group_by(name, type) %>% nest() %>%
    mutate(plt = map(data, \(x) {
      # Colour specification for points
      aes_bee <- if (unique(type) != "prot") {aes(colour = seroclass)} else {NULL}
      # Skip legend if only points of one kind
      param_bee <- if (length(unique(x$seroclass)) == 1) {list(show.legend = F)} else {list(show.legend = T)}
      
      ggplot(x, aes(x = 0, y = scale(value))) %>%
        geom_beebox(aes_bee = aes_bee, param_bee = param_bee, param_box = list(linewidth = 0.5)) +
        labs(x = name, colour = "Serostatus",
             y = ifelse(unique(type) == "prot",
                        "Normalised protein level [AU]",
                        "Relative antibody titer [AU]")) +
        scale_colour_manual(values = c("Positive" = "#1E90FF", "Negative" = "grey")) +
        theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    })) %>%
    ungroup()
  
  # Save in separate files
  file_serol <- plt %>% filter(type != "prot") %>% pull(plt) %>%
    ggsave_multi(file_path = make_fn("data_distribution", "data_distrib_serol.pdf"),
                 n_per_page = 16, return_path = T, return_plots = F, return_blanks = F,
                 nrow = 4, ncol = 4, height = 8, width = 10,
                 guides = "collect", axis_titles = "collect")
  file_prot <- plt %>% filter(type == "prot") %>% pull(plt) %>%
    ggsave_multi(file_path = make_fn("data_distribution", "data_distrib_prot.pdf"),
                 n_per_page = 16, return_path = T, return_plots = F, return_blanks = F,
                 nrow = 4, ncol = 4, height = 8, width = 8,
                 guides = "collect", axis_titles = "collect")
  
  return(c(file_serol, file_prot))
}

# Make a big supplementary table with summary statistics
# - Associations:
#  - Abs vs traits
#  - Proteins vs traits
#  - Proteins vs Ab levels
#  - Proteins vs seroclusters
# - Protein correlation and detectability between methods
# - ROC curve data
# - WGCNA module proteins, ranked by importance 
save_big_table <- function(prot_clin_assoc_tbl, prot_ab_assoc, ifn_prot_assoc, prot_seroclust_assoc,
                           prot_corr, prot_lasso_prop, prot_lasso_others, serol_lasso,
                           prot_module_contrib, wgcna_ukb_comp) {
  # Load protein vs clinical association table
  prot_clin <- read_csv(prot_clin_assoc_tbl, show_col_types = F)
  
  # Get uniprot IDs to add to tables without
  uniprot_id <- prot_clin %>% select(unique_protein_name = protein, uniprot) %>% distinct()
  
  # Combine the lasso variable importance tables into one
  lasso_imp_prot_seroclust <- prot_lasso_prop$varimp %>%
    mutate(class = paste0("serocluster: ", class))
  lasso_imp_prot_others <- imap(prot_lasso_others, \(x, i) {
    x_out <- x$varimp
    # Add class column if not already included, if included add name of variable
    if ("class" %in% colnames(x_out)) {x_out <- mutate(x_out, class = paste0(i, ": ", class))}
    else {x_out <- mutate(x_out, class = i) %>% relocate(class)}
  }) %>% bind_rows()
  lasso_imp_serol <- imap(serol_lasso[1:6], \(x, i) {
    x_out <- x$varimp
    # Add class column if not already included, if included add name of variable
    if ("class" %in% colnames(x_out)) {x_out <- mutate(x_out, class = paste0(i, ": ", class))}
    else {x_out <- mutate(x_out, class = i) %>% relocate(class)}
  }) %>% bind_rows() %>%
    # Change names of serology analytes to match the rest of the tables
    mutate(variable = case_when(str_detect(variable, "Nucleocapsid\\.C") ~ "Anti-Nc",
                                str_detect(variable, "Nucleocapsid_Acro") ~ "Anti-Na",
                                str_detect(variable, "Spike\\.S1S2\\.foldon") ~ "Anti-S1S2",
                                str_detect(variable, "Spike\\.S1") ~ "Anti-S1",
                                str_detect(variable, "RBD") ~ "Anti-RBD",
                                T ~ variable))
  lasso_imp_prot <- bind_rows(lasso_imp_prot_seroclust, lasso_imp_prot_others)
  
  # List of the tables to generate the output file from
  tbl_list <- list(
    # Lasso variable importance (serology)
    "Table S2" = lasso_imp_serol,
    # Protein correlation between technologies
    "Table S3" = prot_corr$corr %>%
      select(protein = name, corr = corr, detectability_alamar = mf_a,
             detectability_olink = mf_o, detected = det_which) %>%
      mutate(detectability_alamar = 1 - detectability_alamar,
             detectability_olink = 1 - detectability_olink) %>%
      rename(unique_protein_name = protein) %>%
      left_join(uniprot_id, by = "unique_protein_name"),
    # Protein antibody associations
    "Table S4" = bind_rows(
      prot_ab_assoc$assoc %>% rename(unique_protein_name = protein),
      ifn_prot_assoc %>% select(unique_protein_name = protein, ab = ifn, estimate, p_value, fdr) %>%
        mutate(ab = paste0("Anti-", ab))
    ) %>% left_join(uniprot_id, by = "unique_protein_name") %>%
      relocate(unique_protein_name, uniprot) %>% arrange(unique_protein_name),
    # Protein serocluster associations (Kruskal-Wallis)
    "Table S5" = prot_seroclust_assoc$kruskal_wallis %>%
      select(-assay_method) %>% rename(unique_protein_name = protein) %>%
      left_join(uniprot_id, by = "unique_protein_name"),
    # Protein serocluster associations (Wilcoxon)
    "Table S6" = prot_seroclust_assoc$wilcoxon %>%
      select(-assay_method) %>% relocate(med_diff, .after = fdr) %>%
      rename(unique_protein_name = protein) %>%
      left_join(uniprot_id, by = "unique_protein_name"),
    # Protein questionnaire associations
    "Table S7" = prot_clin %>% rename(unique_protein_name = protein) %>%
      mutate(protein = str_remove(unique_protein_name, "_\\d$"), .after = unique_protein_name),
    # Lasso variable importance (proteomics)
    "Table S8" = lasso_imp_prot %>%
      left_join(uniprot_id, by = c("variable" = "unique_protein_name")),
    # WGCNA module proteins and importance
    "Table S9" = imap(prot_module_contrib, \(x, i) {
      x$prot_imp %>% select(-PC1) %>% mutate(module = i)
    }) %>% bind_rows(),
    # WGCNA modules and comparing UKB associations
    "Table S10" = wgcna_ukb_comp$test_results %>% select(-trait_num) %>% relocate(module, .after = trait)
  )
  
  # Make a table with column descriptions of the other tables
  tbl_descr <- imap(tbl_list, \(x, i) {
    tibble(sheet = i, column = c(NA, colnames(x)))
  }) %>% bind_rows() %>%
    mutate(description = case_when(column == "unique_protein_name" ~ "Protein assay name where proteins measured repeated times are marked by _1 and _2.",
                                   column == "protein" ~ "Protein assay name.",
                                   column == "clinvar" ~ "Name of questionnaire variable.",
                                   column == "estimate" ~ "Estimate received from regression.",
                                   column == "p_value" ~ "Nominal P-value from statistical test.",
                                   column == "fdr" ~ "P-value adjusted for multiple testing (Benjamini-Hochberg FDR method).",
                                   column == "uniprot" ~ "UniProt identifier for the assayed protein.",
                                   column == "ukb_logp" ~ "Log of P-value for association in UKB. May have multiple if protein was assayed multiple times in the UKB study.",
                                   column == "assoc_serocluster" ~ "Marked with 1 if the protein was associated with serocluster.",
                                   column == "concordant_assoc_abs" ~ "Marked with 1 if the protein was measured multiple times and had concordant associations with an anti-SARS-CoV-2 Ab or anti-IFN AAb.",
                                   column == "ab" ~ "Name of antibody being measured (against SARS-CoV-2 proteins or autoantibody against interferons).",
                                   str_detect(column, "grp\\d") ~ "Groups (seroclusters) being compared.",
                                   column == "med_diff" ~ "Difference in median between the two compared groups (grp1 - grp2).",
                                   column == "corr" ~ "Spearman correlation between the two assays measuring the given protein.",
                                   str_detect(column, "detectability") ~ "Detectability in the different platforms (proportion of samples where the protein was measured above LOD).",
                                   column == "detected" ~ "Which platforms the target was detected in.",
                                   column == "class" ~ "Class being predicted.",
                                   column == "variable" ~ "Name of protein assay (like unique_protein_name).",
                                   column == "importance" ~ "Importance of variable (absolute value of coefficient in Lasso regression).",
                                   column == "sign" ~ "Sign of coefficient in Lasso regression.",
                                   column == "importance_scaled" ~ "Importance of variable scaled between 0 and 1.",
                                   column == "importance_rank" ~ "Rank of importance for WGCNA module eigengene.",
                                   column == "module" ~ "WGCNA module membership.",
                                   column == "trait" ~ "UKB trait.",
                                   column == "n_module" ~ "Number of proteins in module.",
                                   column == "n_other" ~ "Number of proteins in the other modules.",
                                   column == "pos_prop_module" ~ "Proportion of proteins in module that are listed as significantly associated with trait at a Bonferroni-corrected level.",
                                   column == "pos_prop_other" ~ "Proportion of proteins in combined other modules that are listed as significantly associated with trait at a Bonferroni-corrected level.",
                                   column == "lower_ci" ~ "95% confidence interval, lower limit.",
                                   column == "upper_ci" ~ "95% confidence interval, upper limit.",
                                   is.na(column) & sheet == "Table S2" ~ "Variable importance from Lasso regression using antibody levels to predict different traits.",
                                   is.na(column) & sheet == "Table S3" ~ "Protein correlation and detectability in different platforms.",
                                   is.na(column) & sheet == "Table S4" ~ "Protein-antibody associations (linear (logistic) regression).",
                                   is.na(column) & sheet == "Table S5" ~ "Protein-serocluster associations (Kruskal-Wallis test).",
                                   is.na(column) & sheet == "Table S6" ~ "Protein-serocluster associations (Pairwise Wilcoxon rank-sum test).",
                                   is.na(column) & sheet == "Table S7" ~ "Protein-questionnaire associations (linear regression).",
                                   is.na(column) & sheet == "Table S8" ~ "Variable importance from Lasso regression using proteins or antibody levels to predict different traits.",
                                   is.na(column) & sheet == "Table S9" ~ "Protein WGCNA module memberships and their contributions to the module eigengenes (rank of importance).",
                                   is.na(column) & sheet == "Table S10" ~ "WGCNA module associations with UK Biobank (UKB) traits (Fisher exact test)."))
  
  tbl_list <- append(tbl_list, list(description = tbl_descr), after = 0)
  
  fn <- make_fn("master_tables", "results_table.xlsx")
  write_xlsx(tbl_list, fn)
  return(fn)
}
