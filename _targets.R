### Covid vaccination serology pipeline

library(targets)
library(tarchetypes)

pkg <- c("tidyverse", "patchwork", "magrittr", "tidymodels", "vip", "RColorBrewer",
         # Extra options for plotting with ggplot
         "ggh4x", "ggnewscale", "ggbeeswarm", "ggrepel", "ggsignif", "GGally", "ggbiplot", "ggvenn",
         "fpc",                        # For cluster stability function clusterboot()
         "uwot",                       # UMAP
         "ComplexHeatmap", "circlize", # Heatmap visualisation
         "DT",                         # Nice tables for HTML reports
         "table1",                     # Nice tables for static reports
         "knitr",
         "readxl",
         "writexl",
         "ProtPQN",                    # For normalisation with ProtPQN
         "WGCNA",                      # Running WGCNA
         "cluster",                    # Gap statistic with clusGap
         )
tar_option_set(packages = pkg)

# for testing
# Loading packages
# for(p in pkg){library(p, character.only = T)}
# Setting current working directory to where the targets script is
# setwd("scripts/")

# Functions and parameters for pipeline
source("target_params.R")                  # Parameters fed to some targets
source("functions_data.R")                 # Functions for data wrangling targets
source("functions_serology.R")             # Functions for serology targets
source("functions_proteomics.R")           # Functions for proteomics targets
source("functions_interferons.R")          # Functions for interferon targets
source("misc_functions.R")                 # Other small functions called by pipeline steps
source("comp_groups.R")                    # Function for group comparisons
source("tidy_classification.R")            # Function for running Lasso regression

### --- Pipeline steps --- ###

targets_data <- list(
  
  # ## -- Data reading, normalisation -- ##
  # # Read raw serology data
  # tar_target(
  #   serol_files,
  #   c("../data/serol_raw_mfi.csv", "../data/serol_raw_sinfo.csv", "../data/serol_raw_binfo.csv"),
  #   format = "file"
  # ),
  # tar_target(
  #   serol_raw,
  #   read_serol_raw(serol_files, low_conc_s)
  # ),
  # 
  # # Serology normalisation Julia script
  # tar_target(
  #   serol_norm_script,
  #   "ben_mixmod.jl",
  #   format = "file"
  # ),
  # 
  # Normalise serology data for bare-bead correlation
  # (published serology data is this normalised data that has then been scaled)
  # (serol_bbnorm is a list with the elements "mfi": serology_data.csv, "sinfo": sample_info.csv, "binfo": antigen_info.csv)
  tar_target(
    serol_bbnorm,
    #norm_serol(serol_raw, serol_norm_script)
    list(mfi = read_csv("../data/serology_data.csv"),
         sinfo = read_csv("../data/sample_info.csv"),
         binfo = read_csv("../data/antigen_info.csv"))
  ),
  
  # # Plot normalised data
  # tar_target(
  #   serol_norm_plts,
  #   plot_serol_norm(serol_bbnorm),
  #   format = "file"
  # ),
  # 
  # # Read Olink Explore 384 Inflammation panel data for all samples and perform bridge normalisation on two batches
  # tar_target(
  #   olink_binfo,
  #   "../data/olink_binfo.csv",
  #   format = "file"
  # ),
  # tar_target(
  #   olink_batch1,
  #   c("../data/olink_explore_1500_sinfo.csv",
  #     "../data/olink_explore_1500_long_data.csv"),
  #   format = "file"
  # ),
  # tar_target(
  #   olink_batch2,
  #   c("../data/olink_explore_inf_sinfo.csv",
  #     "../data/olink_explore_inf_long_data.csv"),
  #   format = "file"
  # ),
  # tar_target(
  #   olink_bridge_norm,
  #   bridge_normalisation_olink(olink_batch1, olink_batch2, olink_binfo)
  # ),
  # 
  # # Get raw olink data
  # tar_target(
  #   olink_raw,
  #   get_olink_raw(olink_bridge_norm, serol_bbnorm)
  # ),
  # 
  # # Filter Olink data for well-detected proteins
  # tar_target(
  #   olink_filtered,
  #   filter_olink(olink_raw)
  # ),
  # 
  # # Apply ProtPQN normalisation
  # tar_target(
  #   olink_protpqn,
  #   apply_norm_protpqn(olink_filtered)
  # ),
  # 
  # # Detect outliers
  # tar_target(
  #   olink_outlier_det,
  #   olink_protpqn$npm %>%
  #     pivot_longer(cols = -unique_sample_name, names_to = "protein", values_to = "npm") %>%
  #     mutate(panel = "Inflammation") %>%
  #     med_iqr_plt(lab_outl = F)
  # ),
  # 
  # # Data without outliers
  # tar_target(
  #   olink_outl_excl,
  #   outlier_exclude(olink_protpqn, olink_outlier_det)
  # ),
  # 
  # # Read Alamar 200 panel data for all samples
  # tar_target(
  #   alamar_files,
  #   c("../data/alamar_raw_npq.csv", "../data/alamar_raw_sinfo.csv",
  #     "../data/alamar_raw_binfo.csv", "../data/alamar_raw_lod.csv"),
  #   format = "file"
  # ),
  # tar_target(
  #   alamar_raw,
  #   read_alamar_raw(alamar_files)
  # ),
  # 
  # # Filter out proteins with low detectability
  # tar_target(
  #   alamar_filtered,
  #   filter_alamar(alamar_raw)
  # ),
  # 
  # # Apply ProtPQN
  # tar_target(
  #   alamar_protpqn,
  #   apply_norm_protpqn(alamar_filtered)
  # ),
  # 
  # # Detect outliers
  # tar_target(
  #   alamar_outlier_det,
  #   alamar_protpqn$npm %>%
  #     pivot_longer(cols = -unique_sample_name, names_to = "protein", values_to = "npm") %>%
  #     mutate(panel = "Inflammation 200-plex") %>%
  #     med_iqr_plt(lab_outl = F)
  # ),
  # 
  # # Data without outliers
  # tar_target(
  #   alamar_outl_excl,
  #   outlier_exclude(alamar_protpqn, alamar_outlier_det)
  # ),
  # 
  # # Combine all proteomics data into one structure
  # tar_target(
  #   prot_combined,
  #   combine_proteomics(olink_outl_excl, alamar_outl_excl)
  # ),
  
  # Adjust proteomics data for sex
  # (published proteomics data is this normalised and adjusted data, then scaled)
  # (prot_dat is a list with the elements "npm": protein_data.csv, "sinfo": sample_info.csv, "binfo": protein_info.csv)
  tar_target(
    prot_dat,
    #adjust_proteomics(prot_combined)
    list(npm = read_csv("../data/protein_data.csv"),
         sinfo = read_csv("../data/sample_info.csv"),
         binfo = read_csv("../data/protein_info.csv"))
  ),
  
  # Make box + beeswarm plots for all analytes
  tar_target(
    boxplt_all,
    plot_distrib_all(serol_bbnorm, prot_dat, serol_seropositivity, ifn_seropositivity),
    format = "file"
  ),
  
  # Make big table of many summary statistics
  tar_target(
    big_table_file,
    save_big_table(prot_clin_assoc_tbl, prot_ab_assoc, ifn_prot_assoc, prot_seroclust_assoc, prot_corr,
                   prot_lasso_prop, prot_lasso_others, serol_lasso, prot_module_contrib, wgcna_ukb_comp)
  )
)

targets_cohort_summary <- list(
  ## -- Cohort description -- ##
  tar_target(
    demographics_tables,
    make_demographics_tbl(serol_bbnorm)
  ),
  
  tar_target(
    demographics_barplots,
    make_demographics_bar(serol_bbnorm)
  ),
  
  tar_target(
    timing_barplots,
    make_timings_bar(serol_bbnorm)
  ),
  
  tar_render(
    cohort_demographics_report,
    "rmd/cohort_demographics_report.Rmd",
    output_file = make_fn("cohort_demographics", "cohort_demographics_report.html")
  ),
  
  # Age and sex vs seropositivity of N/S/IFN Abs and infection/vaccination
  tar_target(
    age_sex_plots,
    check_age_sex_pos(serol_bbnorm, serol_seropositivity, ifn_seropositivity)
  ),
  
  # Save age/sex plots
  tar_target(
    age_sex_plot_files,
    {
      seropos_fn <- make_fn("cohort_demographics", "age_sex_seropos.pdf")
      vaccinf_fn <- make_fn("cohort_demographics", "age_sex_vaccinf.pdf")
      
      ggsave(seropos_fn, age_sex_plots$plt_abs, height = 12, width = 14)
      ggsave(vaccinf_fn, age_sex_plots$plt_vaccinf, height = 4, width = 6)
      
      c(seropos_fn, vaccinf_fn)
    },
    format = "file"
  )
  
)

targets_serology_analysis <- list(
  ## -- Serology data analysis -- ##
  
  # Boxplots of serology data
  tar_target(
    serol_boxplot,
    make_serol_boxplots(serol_bbnorm)
  ),
  
  # Save PDF version of boxplots
  tar_target(
    serol_boxplot_pdf,
    ggsave_multi(serol_boxplot,
                 file_path = make_fn("multianalyte_serology", "serology_boxplots.pdf"),
                 return_path = T, return_plots = F, n_per_page = 8,
                 nrow = 2, height = 6, width = 13, guides = "collect"),
    format = "file"
  ),
  
  # Seropositivity classification
  tar_target(
    serol_seropositivity,
    classify_serology(serol_bbnorm, n_sd = 6)
  ),
  
  # Summary tables
  tar_target(
    serol_table,
    make_serol_table(serol_bbnorm, serol_seropositivity)
  ),
  
  # Serology report
  tar_render(
    multianalyte_serology_report,
    "rmd/multianalyte_serology_report.Rmd",
    output_file = make_fn("multianalyte_serology", "multianalyte_serology_report.html")
  ),
  
  # Hierarchical clustering of serology data
  tar_target(
    seroclust,
    run_serol_hclust(serol_bbnorm)
  ),
  
  # Table over serocluster characteristics
  tar_target(
    seroclust_table,
    make_seroclust_table(seroclust, serol_bbnorm, serol_seropositivity, ifn_seropositivity)
  ),
  
  # Characterising "mismatching" samples
  tar_target(
    seroclust_mismatch,
    get_seroclust_mismatch(seroclust, serol_bbnorm)
  ),
  
  # Heatmap of clustering
  tar_target(
    seroclust_hm,
    make_seroclust_heatmap(seroclust, clust_cols)
  ),
  
  # Boxplots of Abs in clusters and questionnaire annotation
  tar_target(
    seroclust_boxplt,
    make_seroclust_boxplt(seroclust)
  ),
  
  # Response over time plots in clusters
  tar_target(
    seroclust_response,
    make_seroclust_responseplt(seroclust, serol_bbnorm, ab_cols)
  ),
  
  tar_render(
    seroclust_report,
    "rmd/serocluster_comparison_report.Rmd",
    output_file = make_fn("serocluster_comparison", "serocluster_comp_report.html")
  ),
  
  # Figure summarising serology analyses
  tar_target(
    serol_figure,
    make_serology_figure(serol_boxplot, seroclust_hm, seroclust_mismatch, seroclust_response),
    format = "file"
  ),
  
  # Classification using serology
  tar_target(
    serol_lasso,
    serol_classification(serol_bbnorm, seroclust, use_ifn = F, clust_cols2)
  ),
  # Save figure
  tar_target(
    serol_lasso_fig,
    {
      imap(serol_lasso[c("age_group", "sex", "region", "serocluster", "vaccinated", "infected")], \(x, i) {
        x$roc_curve +
          labs(title = paste0(str_to_sentence(str_replace(i, "_", " ")), ", AUC: ",
                              x$metrics %>%
                                filter(.metric == "roc_auc") %>%
                                pull(.estimate) %>% signif(2))) +
          theme_classic(12)
      }) %>% wrap_plots(axis_titles = "collect") %>%
        ggsave(filename = make_fn("multianalyte_serology", "serol_lasso.pdf"),
               height = 7, width = 12)
    },
    format = "file"
  )
  
)

targets_proteomics_analysis <- list(
  ## -- Proteomics data analysis -- ##
  
  # # Protein correlations between technologies
  # tar_target(
  #   prot_corr,
  #   run_prot_corr(olink_raw, alamar_raw, olink_outlier_det, alamar_outlier_det)
  # ),
  # 
  # # Report for protein correlation plots
  # tar_render(
  #   prot_corr_report,
  #   "rmd/prot_corr_report.Rmd",
  #   output_file = make_fn("protein_correlation", "prot_corr_report.html")
  # ),
  # 
  # # Make some QC plots
  # # Protein detectability
  # tar_target(
  #   prot_detectability,
  #   make_prot_detectability_plot(olink_raw, alamar_raw)
  # ),
  # 
  # # Protein CV
  # tar_target(
  #   prot_cv,
  #   make_prot_cv_plot(olink_raw, alamar_raw, prot_combined)
  # ),
  
  # Protein IQR
  tar_target(
    prot_iqr,
    make_prot_iqr_plot(prot_dat)
  ),
  
  # Sample PCA
  tar_target(
    prot_pca,
    make_prot_pca_plot(prot_dat, clust_cols)
  ),
  # Sample UMAP
  tar_target(
    prot_umap,
    make_prot_umap_plot(prot_dat)
  ),
  
  # # Figure of QC plots + outliers
  # tar_target(
  #   prot_qc_plot,
  #   make_prot_qc_plot(outliers_a = alamar_outlier_det, outliers_o = olink_outlier_det,
  #                     prot_det = prot_detectability, prot_iqr = prot_iqr, prot_cv = prot_cv, prot_pca = prot_pca)
  # ),
  # 
  # # Save protein QC plot
  # tar_target(
  #   prot_qc_file,
  #   ggsave(filename = make_fn("protein_qc", "protein_qc_plot.pdf"), plot = prot_qc_plot,
  #          height = 9, width = 10)
  # ),
  
  # Associations of proteins with questionnaire variables
  tar_target(
    prot_clin_assoc,
    run_prot_clin_assoc(prot_combined, ukb_assoc1)
  ),
  
  # Save a table of the association results
  tar_target(
    prot_clin_assoc_tbl,
    save_prot_clin_assoc(prot_clin_assoc, prot_seroclust_assoc, ifn_prot_hms),
    format = "file"
  ),
  
  # Associations for samples in serocluster 1
  tar_target(
    prot_clin_assoc_sc1,
    {
      # Prepare data with only samples from serocluster 1
      in_temp <- prot_combined
      in_temp$sinfo <- in_temp$sinfo %>%
        left_join(seroclust$serocluster %>% select(unique_sample_name, serocluster),
                  by = "unique_sample_name") %>%
        filter(serocluster == 1)
      in_temp$npm <- in_temp$npm %>% filter(unique_sample_name %in% in_temp$sinfo$unique_sample_name)
      # Run analysis
      run_prot_clin_assoc(in_temp, ukb_assoc1)
    }
  ),
  # Save table
  tar_target(
    prot_clin_assoc_sc1_tbl,
    save_prot_clin_assoc(prot_clin_assoc_sc1, prot_seroclust_assoc, ifn_prot_hms, "prot_clin_assoc_sc1.csv")
  ),
  # Plot results
  tar_target(
    prot_clin_assoc_sc1_plt,
    ggsave(filename = make_fn("protein_clinical_association", "prot_clin_assoc_sc1.pdf"),
           plot = plot_prot_clin_assoc(prot_clin_assoc_sc1_tbl), height = 5, width = 10)
  ),
  
  # Protein-questionnaire associations per region
  tar_target(
    prot_clin_assoc_region,
    sapply(c("Stockholm", "Gothenburg"), \(x) {
      prot_in <- prot_combined
      prot_in$sinfo <- prot_in$sinfo %>% filter(region == x)
      prot_in$npm <- prot_in$npm %>% filter(unique_sample_name %in% prot_in$npm$unique_sample_name)

      run_prot_clin_assoc(prot_in, ukb_assoc1)
    }, simplify = F, USE.NAMES = T)
  ),
  
  # Plot the betas per region
  tar_target(
    prot_clin_beta_comparison,
    make_clin_beta_plot(prot_clin_assoc_region)
  ),
  
  # Associations of Abs with proteins
  tar_target(
    prot_ab_assoc,
    run_prot_ab_assoc(prot_dat, serol_seropositivity)
  ),
  
  # Save heatmaps of associations (full size and selected proteins)
  tar_target(
    prot_ab_hm,
    {
      ggsave(make_fn("protein_ab_association", "prot_ab_assoc_heatmap.pdf"),
             (prot_ab_assoc$heatmap_fdr + labs(tag = "A")) +
               (prot_ab_assoc$heatmap_p + labs(tag = "B")) +
               plot_layout(guides = "collect"),
             height = 12, width = 17)
    },
    format = "file"
  ),
  
  # Report for protein-ab associations
  tar_render(
    prot_ab_assoc_report,
    "rmd/prot_ab_assoc_report.Rmd",
    output_file = make_fn("protein_ab_association", "prot_ab_assoc_report.html")
  ),
  
  # Assocation of clusters with proteins
  tar_target(
    prot_seroclust_assoc,
    run_prot_seroclust_assoc(prot_dat, seroclust)
  ),
  tar_target(
    prot_seroclust_plts,
    make_prot_seroclust_plots(prot_seroclust_assoc, prot_dat, seroclust, clust_cols2)
  ),
  
  # Save boxplots
  tar_target(
    prot_seroclust_boxplts,
    save_prot_seroclust_boxplts(prot_seroclust_plts)
  ),
  
  # Lasso analysis to predict serocluster with a 70:30 train:test split
  tar_target(
    prot_lasso_prop,
    run_prot_ml(prot_dat, seroclust, .7, clust_cols2)
  ),
  
  # Lasso for predicting sex, age, and region with proteins (before adjusting for sex)
  tar_target(
    prot_lasso_others,
    run_prot_ml_others(prot_combined, .7)
  ),
  
  # Lasso for predicting vaccination on age-adjusted protein levels
  tar_target(
    prot_lasso_vacc,
    run_prot_ml_vacc(prot_combined, .7)
  ),
  
  # Put together a figure with the results from the region-split lasso
  tar_target(
    prot_lasso_fig,
    {
      ggsave(plot = wrap_plots(append(list(prot_lasso_prop$roc_curve +
                                             labs(title = "a", colour = "Serocluster") +
                                             annotate("text", x = 0.75, y = 0.1, size = 6,
                                                      label = prot_lasso_prop$metrics %>%
                                                        filter(.metric == "roc_auc") %>%
                                                        pull(.estimate) %>% signif(2) %>% paste0("AUC: ", .))),
                                      imap(prot_lasso_others[c("age_group", "sex", "region", "vaccinated", "infected")] %>% unname(), \(x, i) {
                                        x$roc_curve + labs(title = letters[i + 1],
                                                           colour = switch(i, "1" = "Age group")) +
                                          annotate("text", x = 0.75, y = 0.1, size = 6,
                                                   label = x$metrics %>% filter(.metric == "roc_auc") %>%
                                                     pull(.estimate) %>% signif(2) %>% paste0("AUC: ", .))
                                      })), nrow = 3, axis_titles = "collect") &
               coord_cartesian(clip = "off") &
               theme(plot.title = element_text(hjust = 0)),
             filename = make_fn("protein_lasso", "prot_lasso_fig.pdf"),
             height = 9, width = 10)
    },
    format = "file"
  ),
  
  tar_render(
    prot_lasso_report,
    "rmd/prot_lasso_report.Rmd",
    output_file = make_fn("protein_lasso", "prot_lasso_report.html")
  ),
  
  # Look at CRP levels
  tar_target(
    crp_positivity,
    classify_crp(prot_dat, n_sd = c(3, 6)) %>%
      mutate(class = case_when(class == "Negative" ~ "Low",
                               class == ">3 SD" ~ "Medium",
                               class == ">6 SD" ~ "High")) %>%
      mutate(class = factor(class, levels = c("Low", "Medium", "High")))
  ),
  
  # Plot CRP levels
  tar_target(
    crp_plot,
    make_crp_plot(crp_positivity, serol_bbnorm, palette = "Blues"),
    format = "file"
  ),
  
  # Characterise CRP groups
  tar_target(
    crp_char,
    characterise_crp(crp_positivity, prot_dat, seroclust, prot_clust),
    format = "file"
  ),
  
  # WGCNA
  tar_target(
    prot_wgcna,
    run_wgcna(prot_dat)
  ),
  
  # Visualising WGCNA results
  tar_target(
    wgcna_vis,
    visualise_wgcna_modules(prot_dat, prot_wgcna, seroclust)
  ),
  
  # Save table of WGCNA modules and their proteins
  tar_target(
    wgcna_module_proteins,
    {
      tibble(protein = names(prot_wgcna$wgcna$colors), module = prot_wgcna$wgcna$colors) %>%
        left_join(prot_dat$binfo %>% select(unique_protein_name, assay, uniprot),
                  by = c("protein" = "unique_protein_name")) %>%
        relocate(assay) %>% rename(unique_assay_name = protein) %>%
        arrange(module)
    }
  ),
  tar_target(
    wgcna_module_table,
    {
      fn <- make_fn("protein_clustering", "module_proteins.csv")
      write_csv(x = wgcna_module_proteins, file = fn)
      fn
    },
    format = "file"
  ),
  
  # Contributions of proteins to the MEs of each WGCNA module
  tar_target(
    prot_module_contrib,
    wgcna_module_pca(prot_dat, prot_wgcna, prot_clust_effect)
  ),
  
  # Make plot of how proteins end up in the WGCNA modules, also looking at overlapping proteins if measured with both platforms
  tar_target(
    wgcna_prot_overlap,
    check_wgcna_overlap(wgcna_module_proteins)
  ),
  
  # Save the plot
  tar_target(
    wgcna_prot_overlap_file,
    ggsave(make_fn("protein_clustering", "module_prot_assignment.pdf"),
           wgcna_prot_overlap$plt, height = 10, width = 7)
  ),
  
  # Load UK Biobank associations data
  # Supplementary tables 5 and 7 of Sun et al. https://doi.org/10.1038/s41586-023-06592-6,
  # loaded into R to have columns like age_beta, age_se, age_logp, sex_beta, etc
  # as well as the 'UKBPPP ProteinID' column split into 'assay', 'uniprot', 'olinkid', and 'version'.
  # 'Protein name' is named 'protein_long'
  # Associations with age, sex, BMI
  tar_target(
    ukb_assoc1,
    read_csv("../data/ukb_st5.csv")
  ),
  # Associations with medication and disease
  tar_target(
    ukb_assoc2,
    read_csv("../data/ukb_st7.csv")
  ),
  
  # Associations per module
  tar_target(
    wgcna_ukb_assoc,
    make_wgcna_ukb_table(prot_dat, prot_wgcna, ukb_assoc1, ukb_assoc2)
  ),
  
  # Save a plot of the associations
  tar_target(
    wgcna_ukb_plot,
    {
      plt <- make_wgcna_ukb_plots(wgcna_ukb_assoc)
      ggsave(make_fn("protein_clustering", "wgcna_ukb_assoc.pdf"), plt, height = 6, width = 10)
    },
    format = "file"
  ),
  
  # Make a heatmap of the WGCNA UKB assocations
  tar_target(
    wgcna_ukb_comp,
    make_wgcna_ukb_heatmap(wgcna_ukb_assoc)
  ),
  # Save heatmap
  tar_target(
    wgcna_ukb_heatmap,
    ggsave(make_fn("protein_clustering", "wgcna_ukb_heatmap.pdf"), wgcna_ukb_comp$plot, height = 14, width = 12)
  ),
  
  # Load HPA annotations for proteins
  tar_target(
    prot_hpa_annot,
    read_xlsx("../data/protein_hpa_annotation.xlsx", sheet = 1) %>%
      select(assay = Gene, contains("Blood cell")) %>%
      rename_with(\(x) {
        case_when(str_detect(x, "Specificity") ~ "blood_cell_spec",
                  str_detect(x, "Enriched") ~ "blood_cell_enrich",
                  str_detect(x, "Distribution") ~ "blood_cell_distrib",
                  T ~ x)
      })
  ),
  
  # Plot blood cell enrichment in modules
  tar_target(
    prot_module_blood_plt,
    make_wgcna_hpa_plot(wgcna_module_proteins, prot_hpa_annot)
  ),
  tar_target(
    prot_module_blood_plt_file,
    ggsave(filename = make_fn("protein_clustering", "wgcna_blood_enrich.pdf"),
           plot = prot_module_blood_plt, height = 5, width = 8),
    format = "file"
  ),
  
  # Clustering samples based on WGCNA MEs
  tar_target(
    prot_clust,
    run_clustering_me(prot_wgcna, linearise_data = F, scale_data = T, rseed = 123,
                      clust_distance = "manhattan", clust_method = "ward.D2",
                      gap_b = 100, nrep = 50, verbose = F)
  ),
  
  # Plot PC1-4 coloured by proteocluster
  tar_target(
    prot_clust_pca,
    make_prot_pca1to4(prot_dat, prot_clust)
  ),
  # Save the plot
  tar_target(
    prot_clust_pca_file,
    ggsave(filename = make_fn("protein_clustering", "prot_clust_pca.pdf"),
           plot = prot_clust_pca, height = 7, width = 7),
    format = "file"
  ),
  
  # Make PCA of the module eigengenes
  tar_target(
    prot_me_pca,
    make_prot_pca_biplot(prot_wgcna, prot_clust)
  ),
  tar_target(
    prot_me_pca_file,
    ggsave(filename = make_fn("protein_clustering", "prot_me_pca.pdf"),
           plot = prot_me_pca, height = 6, width = 8),
    format = "file"
  ),
  
  # Characterise proteoclusters
  tar_target(
    prot_clust_char,
    characterise_proteoclust(prot_wgcna, prot_clust, prot_dat, seroclust, crp_positivity, ifn_seropositivity)
  ),
  
  # Calculate effect size (difference in mean protein level) between proteotypes
  tar_target(
    prot_clust_effect,
    compute_clust_effect(prot_dat$npm, prot_clust, prot_wgcna, clust_vs_rest = T) %>%
      mutate(module = factor(module, levels = rev(append(setdiff(sort(unique(module)), "grey"), "grey"))))
  ),
  
  # Calculate mean effect sizes between proteotypes on the MEs themselves instead of constituent proteins
  tar_target(
    prot_me_effect,
    compute_clust_effect(prot_wgcna$wgcna$MEs %>% rownames_to_column("unique_sample_name") %>%
                           rename_with(\(x) str_to_sentence(str_remove(x, "ME")), .cols = contains("ME", ignore.case = F)),
                              prot_clust, prot_wgcna, clust_vs_rest = T) %>%
      mutate(comparison = as.character(comparison),
             name = factor(name, levels = rev(append(setdiff(sort(unique(name)), "Grey"), "Grey"))))
  ),
  tar_target(
    prot_me_effect_plts,
    make_prot_clust_effect_plots(prot_me_effect, spread = 0.6, point_size = 0.5)
  ),
  tar_target(
    prot_me_effect_file,
    ggsave(make_fn("protein_clustering", "prot_effect_me.pdf"), prot_me_effect_plts, width = 6, height = 5.5),
    format = "file"
  ),
  
  tar_target(
    prot_clust_fig,
    save_proteotype_fig(prot_me_pca, prot_clust_char, prot_me_effect_plts,
                        prot_clust_effect, prot_seroclust_plts, wgcna_prot_overlap),
    format = "file"
  ),
  
  # Report on proteoclusters/proteotypes
  tar_render(
    prot_clust_report,
    "rmd/prot_clustering_report.Rmd",
    output_file = make_fn("protein_clustering", "prot_clustering_report.html")
  )
  
)

targets_interferon_analysis <- list(
  ### --- Interferon analysis --- ###
  
  # Get interferon data
  tar_target(
    ifn_dat,
    make_ifn_data(serol_bbnorm)
  ),
  
  # Classify serology of interferons
  tar_target(
    ifn_seropositivity,
    classify_ifn(ifn_dat, n_sd = 12)
  ),
  
  # Make boxplots per cluster
  tar_target(
    ifn_vs_seroclust,
    make_ifn_clust_plots(ifn_dat, ifn_seropositivity, seroclust)
  ),
  
  # Make summary table
  tar_target(
    ifn_table,
    make_ifn_table(ifn_dat, ifn_seropositivity, seroclust)
  ),
  
  # Interferon association with age and sex
  tar_target(
    ifn_assoc,
    check_ifn_assoc(ifn_dat, ifn_seropositivity, serol_seropositivity)
  ),
  
  # Make summary plot and table for IFN associations
  tar_target(
    ifn_assoc_summary,
    make_ifn_assoc_summary(ifn_assoc)
  ),
  
  # Save interferon summary
  tar_target(
    ifn_assoc_summary_plot,
    save_ifn_summary(ifn_assoc_summary),
    format = "file"
  ),
  
  # Association between IFN positivity and protein levels
  tar_target(
    ifn_prot_assoc,
    check_ifn_prot_assoc(ifn_seropositivity, prot_dat)
  ),
  
  # Make heatmaps of the results
  tar_target(
    ifn_prot_hms,
    make_ifn_prot_hms(ifn_prot_assoc, prot_ab_assoc)
  ),
  
  # Save the IFN protein association plots
  tar_target(
    ifn_prot_assoc_p,
    {
      filename <- make_fn("interferon_analysis", "ifn_prot_assoc_p.pdf")
      pdf(file = filename, height = 7, width = 5.5)
      draw(ifn_prot_hms$heatmap_p)
      dev.off()
      filename
    }
  ),
  tar_target(
    ifn_prot_assoc_fdr,
    ggsave(filename = make_fn("interferon_analysis", "ifn_prot_assoc_fdr.pdf"),
           plot = ifn_prot_hms$heatmap_fdr,
           height = 12, width = 8)
  ),
  
  # Make Rmd report
  tar_render(
    interferon_report,
    "rmd/interferon_report.Rmd",
    output_file = make_fn("interferon_analysis", "interferon_report.html")
  )
)

targets_package_versions <- list(
  tar_target(
    pkg_versions,
    {
      # Create dependency on package list
      pkgs <- pkg
      pkg_version_tbl()
    }
  ),
  
  # Save a table
  tar_target(
    pkg_versions_file,
    {
      fn <- make_fn("package_versions", "package_versions.csv")
      write_csv(pkg_versions, fn)
      fn
    },
    format = "file"
  )
)

# Run pipeline
list(targets_data,
     targets_cohort_summary,
     targets_serology_analysis,
     targets_proteomics_analysis,
     targets_interferon_analysis,
     targets_package_versions)
