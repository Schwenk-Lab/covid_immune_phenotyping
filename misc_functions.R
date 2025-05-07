# Miscellaneous functions used by pipeline step (target) functions

# Change theme_classic so that axis text is black instead of dark grey
theme_classic <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"))
}

# Overwrite name clashes from packages such as MASS and many bioconductor packages
select <- dplyr::select
rename <- dplyr::rename
rbind <- base::rbind
intersect <- dplyr::intersect
slice <- dplyr::slice

# Make a filename for saving a file
make_fn <- function(results_dir = "", name = "") {
  paste0(paste0("../results/", results_dir, "/"),
         format(Sys.time(), "%Y-%m-%d_%H%M%S"),
         "_", name)
}

# Convenience function for picking out only the SARS-CoV-2-related binders for serology
get_cov <- function(x) {
  select(x, unique_sample_name, contains(c("anti_n", "anti_s", "rbd")))
}

# Function for renaming serology binders in long format data for plotting
rename_serol <- function(x, anti_all = T) {
  case_when(str_detect(x, "anti_nc") ~ "Anti-Nc",
            str_detect(x, "anti_na") ~ "Anti-Na",
            str_detect(x, "anti_s1s2") ~ "Anti-S1S2",
            str_detect(x, "anti_s1") ~ "Anti-S1",
            str_detect(x, "anti_rbd") ~ "Anti-RBD",
            str_detect(x, "anti_human_igg") ~ "Anti-human IgG",
            str_detect(x, "anti_human_igm") ~ "Anti-human IgM",
            str_detect(x, "anti_human_iga") ~ "Anti-human IgA",
            str_detect(x, "anti_rabbit_igm") ~ "Anti-rabbit IgM",
            x == "unique_sample_name" ~ x,
            T ~ paste0(ifelse(anti_all, "Anti-", ""), x))
}

# Hierarchical clustering function for stability analysis by bootstrapping
boot_clustfun <- function(x, dist_method = clust_dist, cl_method = clust_method, n_cl = n_clust) {
  # x, matrix with data to cluster
  
  clust_results <- dist(x, method = dist_method) %>% hclust(method = cl_method)
  clusts <- cutree(clust_results, n_cl)
  
  out_list <- list("result" = clust_results,
                   "nc" = length(unique(clusts)),
                   "clusterlist" = lapply(sort(unique(clusts)), function(x) {clusts == x}),
                   "partition" = clusts,
                   "clustermethod" = "hclust")
}

# Compute a proportion confidence interval (two-tailed)
prop_ci <- function(p, N, conf_level = 0.95) {
  z_crit <- qnorm((1 - conf_level) / 2, lower.tail = F)
  ci <- c(p - z_crit * sqrt((p * (1 - p)) / N),
          p + z_crit * sqrt((p * (1 - p)) / N))
  # If negative proportion substitute for 0
  ci[ci < 0] <- 0
  
  return(ci)
}

# Density cutoff
density_cutoff <- function(x, n_sd, neg_prop, conf_level = 0.95) {
  # Get population density peak and distance from peak
  dens <- density(x)
  peak <- dens$x[which.max(dens$y)]
  pdist <- abs(x - peak)
  
  # Get variation of expected negative proportion of samples around peak
  qtl <- quantile(pdist, probs = neg_prop)
  variation <- sd(x[pdist <= qtl])
  
  # Classify points
  cutoff <- peak + n_sd * variation
  cls <- x > cutoff
  prop <- sum(cls) / length(cls)
  
  # Confidence interval for proportion
  ci <- prop_ci(prop, length(cls), conf_level = conf_level)
  
  return(list("classes" = cls,
              "cutoff" = cutoff,
              "pos_prop" = prop,
              "confint" = ci,
              "conf_level" = conf_level))
}

# Combine replicates if enough volume
combine_repl <- function(dat_in) {
  npm_merge <- dat_in$npm %>%
    # Add sample id to group by
    left_join(dat_in$sinfo %>% select(unique_sample_name, sample_id, volume), by = "unique_sample_name") %>%
    filter(volume > 10) %>% select(-volume) %>%
    pivot_longer(cols = -c(unique_sample_name, sample_id), names_to = "assay", values_to = "value") %>%
    # Compute mean per sample and protein, merging replicates and leaving other samples intact
    group_by(sample_id, assay) %>%
    summarise(value = mean(value), .groups = "keep") %>%
    pivot_wider(id_cols = sample_id, names_from = assay, values_from = value) %>%
    rename(unique_sample_name = sample_id) %>%
    ungroup()
  
  # Remove duplicated entries in sample info and remove the replicate labels (turning sample_id and unique_sample_name into identical columns)
  sinfo_out <- dat_in$sinfo %>%
    filter(volume > 10) %>%
    filter(!duplicated(sample_id)) %>%
    mutate(unique_sample_name = str_remove(unique_sample_name, "_[:alpha:]$"))
  
  # For consistency, keep same order of samples in both
  npm_merge <- npm_merge[match(sinfo_out$unique_sample_name, npm_merge$unique_sample_name), ]
  
  dat_out <- dat_in
  dat_out$npm <- npm_merge
  dat_out$sinfo <- sinfo_out
  
  return(dat_out)
}

# Calculating CV%
cv <- function(x, na.rm = F) {
  sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm)) * 100
}

# Function for making PCA in different contexts (simple two PC scatter, >2 PC pairs plot, two-PC biplot)
plot_pca <- function(x_in, pcs = c(1, 2), col_vec = NULL, col_scale = NULL, fill_scale = NULL,
                    scale = T, center = T, biplot = F, biplot_ellipse = T) {
  # Run PCA
  pca <- x_in %>%
    prcomp(scale. = scale, center = center) %>%
    # Use `summary()` to get proportion of variance for PCs unless biplot
    # (does not accept summary.prcomp objects)
    (\(x) {if (biplot & length(pcs) == 2) x else summary(x)})()
  
  if (biplot & length(pcs) == 2) {
    # Make axis names with proportions of variance
    axis_names <- summary(pca) %>%
      pluck("importance") %>%
      as.data.frame() %>%
      slice(2) %>%
      select(all_of(pcs)) %>%
      unlist() %>%
      signif(2)
    axis_names <- paste0("PC", pcs, " (", axis_names * 100, "%)")
    
    plt_out <- ggbiplot(pca, choices = pcs, groups = col_vec, ellipse = biplot_ellipse, varname.size = 5) +
      col_scale + fill_scale +
      labs(x = axis_names[1], y = axis_names[2], colour = "Group", fill = "Group") +
      theme_classic(16) +
      theme(aspect.ratio = 1)
    
  } else if (biplot & length(pcs) != 2) {
    stop("The biplot needs exactly two PCs.")
    
  } else {
    # Make data frame for plotting, unless making a biplot
    plt_dat <- pca %>%
      pluck("x") %>%
      as.data.frame() %>%
      select(all_of(pcs)) %>%
      # Add importance to names
      rename_with(\(x) {
        paste0(x, " (", signif(pca$importance[2, pcs], 2) * 100, "%)")
      }) %>%
      (\(x) {if (!is.null(col_vec)) {
        mutate(x, colr = col_vec)
      } else {x}})()
    
    # Get the names (used if 2 PCs) and indices (used if more PCs) of the columns to plot
    cols_to_plot <- str_subset(colnames(plt_dat), paste(paste0("^PC", pcs), collapse = "|"))
    inds_to_plot <- str_which(colnames(plt_dat), paste(paste0("^PC", pcs), collapse = "|"))
    
    # Make plot
    if (length(pcs) == 2) {
      # Two PCs
      plt_out <- ggplot(plt_dat, aes(x = .data[[cols_to_plot[1]]], y = .data[[cols_to_plot[2]]],
                                     colour = if (!is.null(col_vec)) {colr} else {NULL})) +
        geom_point(size = 2.5, alpha = 0.8) +
        col_scale +
        labs(colour = "Group") +
        theme_classic(16)
      
    } else if (length(pcs) > 2) {
      # Pairs plot if more than two PCs
      plt_out <- ggpairs(plt_dat, mapping = aes(colour = if (!is.null(col_vec)) {colr} else {NULL}),
                         columns = inds_to_plot,
                         lower = list(continuous = wrap("points", alpha = 0.5)),
                         diag = list(continuous = wrap("densityDiag", alpha = 0.5))) +
        theme_bw() +
        theme(strip.background = element_rect(fill = "white"),
              panel.grid = element_blank())
      
      if (!is.null(col_scale)) {plt_out <- plt_out + col_scale + labs(colour = "Group")}
      if (!is.null(fill_scale)) {plt_out <- plt_out + fill_scale + labs(fill = "Group")}
      
    } else {
      stop("Select at least two PCs to plot.")
    }
  }
  
  return(plt_out)
}

ggsave_multi <- function(...,
                         file_path = NULL,
                         n_per_page = NULL,
                         n_pages = NULL,
                         return_path = T,
                         return_plots = F,
                         return_blanks = F,
                         height = 7,
                         width = 7,
                         nrow = NULL,
                         ncol = NULL,
                         guides = NULL,
                         axis_titles = NULL) {
  # Function for saving multiple ggplot plots (possibly) on multiple pages into a pdf file.
  # If no file path is provided the function only generates a nested list where the provided plots
  # have been distributed into
  # ..., either a list of ggplots or ggplots separately, to be plotted multiple together
  # file_path, character of length 1 that, if provided, will be used as the file name of the output file (PDF file with multiple plots per page)
  # n_per_page, numeric indicating the number of plots to include on each page
  # n_pages, numeric indicating how many pages the plots should be distributed on
  # May give fewer pages than specified if the specified number of pages would result in more than one page getting fewer plots than the previous pages
  # n_per_page and n_pages are mutually exclusive, provide one of them
  # return_path, logical, if TRUE the file path will be returned. File path is not returned if return_plots is TRUE
  # return_plots, logical, if TRUE the function returns the plots in a nested list, each element containing the plots plotted on each page
  # return_blanks, logical, if TRUE the output will contain plot_spacer() elements in pages with fewer real plots, to fill up the space on the page
  # height, width, numerics of length 1 to be used as input to the pdf function
  # nrow, ncol, numerics given to the wrap_plots function, specifying the number of rows and columns, respectively, to have on each page
  # guides, character given to the wrap_plots function, specifying how legends should be treated
  # axis_titles, character given to wrap_plots function, specifying if axis titles should be bunched together
  
  if (all(sapply(list(...), is.ggplot))) {
    # Plots not given as a list, put in a list
    plts <- list(...)
    
  } else if (length(list(...)) > 1) {
    # More than one unnamed argument but not all are ggplots?! no output!
    stop("As unnamed arguments, provide ggplot objects or one list containing ggplot objects")
    
  } else if (is.list(..1)) {
    # Plots are in a list (first argument of ... only)
    plts <- ..1
    
  } else {
    # No plots, no list, no output!
    stop("As unnamed arguments, provide ggplot objects or one list containing ggplot objects")
    
  }
  
  # If both n_per_page and n_pages are specified, complain
  # Also complain if neither is specified
  if ((!is.null(n_per_page) & !is.null(n_pages)) |
      (is.null(n_per_page) & is.null(n_pages))) {
    stop("Specify either n_per_page or n_pages")
  }
  
  # Split list of plots into list with multiple plots per element
  if (!is.null(n_per_page)) {
    splits <- split(seq_along(plts),
                    ceiling(seq_along(plts) / n_per_page))
    
  } else if (!is.null(n_pages)) {
    splits <- split(seq_along(plts),
                    ceiling(seq_along(plts) / ceiling(length(plts) / n_pages)))
  }
  
  plts_to_plot <- lapply(splits, \(x) {
    # Count to see if the page has fewer plots
    if (length(x) < max(length(splits))) {
      # If so, fill up with plot spacers
      return(c(plts[x], rep(list(plot_spacer()), max(lengths(splits)) - length(x))))
      
    } else {
      return(plts[x])
      
    }
  })
  
  # Save plots as a pdf file if a path has been provided
  if (!is.null(file_path)) {
    pdf(file_path, onefile = T, height = height, width = width)
    for (i in plts_to_plot) {
      print(
        wrap_plots(i, guides = guides, nrow = nrow, ncol = ncol, axis_titles = axis_titles)
      )
    }
    dev.off()
  }
  
  # Return the plot list if desired
  if (return_plots) {
    if (return_blanks) {return(plts_to_plot)}
    else {return(lapply(splits, \(x) {plts[x]}))}
  } else if (return_path) {
    return(file_path)
  } else {return(NULL)}
}

# Gap statistic for optimal number of clusters in hierarchical clustering
# Optimal number of clusters according to gap statistic
get_optimal_k <- function(x, k_max = 15, b = 50, distance = "manhattan", method = "ward.D2", verbose = F) {
  hclust2 <- function(x, k, distance, method) {
    clust <- x %>% dist(method = distance) %>% hclust(method = method) %>% cutree(k)
    return(list("cluster" = clust))
  }
  gap_stat <- clusGap(x, FUNcluster = hclust2, K.max = k_max, B = b,
                      distance = distance, method = method, verbose = verbose)
  
  # Get optimal k as smallest k where gap_k >= gap_k+1 - SE_k+1
  gap_tab <- gap_stat$Tab %>%
    as.data.frame() %>%
    # Gap minus SE for next k
    mutate(gms = lead(gap) - lead(SE.sim),
           gap_larger = gap >= gms)
  
  optimal_k <- which(gap_tab$gap == (gap_tab %>% filter(gap_larger) %>% pull(gap) %>% min()))
  
  return(list("results" = gap_stat, "optimal_k" = optimal_k))
}

#' Wrapper function for combining geom_beeswarm and geom_boxplot in ggplot2
#'
#' @param gg_in ggplot object to which the beeswarm and boxplot layers should be added
#' @param aes_bee `aes()` call for the beeswarm to overwrite defaults or add new arguments
#' @param param_bee Named list containing extra parameters for `geom_beeswarm()`
#' @param aes_box `aes()` call for the boxplot to overwrite defaults or add new arguments
#' @param param_box Named list containing extra parameters for `geom_boxplot()`
#' @param new_col_scale Logical indicating if a new colour scale should be added between the bees and the boxes. Only required if colour is specified in aes() for both the beeswarm and the boxplot.
#' @param new_fill_scale Logical like `new_col_scale` but for fill scales
#' @param col_scale_bee Call to a `scale_colour_x` function to define colour scale for the beeswarm plot
#' @param fill_scale_bee Call to a `scale_fill_x` function to define colour scale for the beeswarm plot
#' @param col_scale_box Like `col_scale_bee` but for the boxplots
#' @param fill_scale_box Like `fill_scale_bee` but for the boxplots
#'
#' @return A ggplot object with added beeswarm and boxplot layers.
#' @export
#'
#' @examples ggplot(iris, aes(x = Species, y = Sepal.Length)) %>%
#'             geom_beebox(aes_bee = aes(colour = Species),
#'             param_bee = list(size = 3),
#'             param_box = list(width = 0.05)) +
#'             scale_colour_brewer(palette = "Pastel1") +
#'             theme_classic()
geom_beebox <- function(gg_in,
                        aes_bee = NULL, param_bee = NULL,
                        aes_box = NULL, param_box = NULL,
                        new_col_scale = F, new_fill_scale = F,
                        col_scale_bee = NULL, fill_scale_bee = NULL,
                        col_scale_box = NULL, fill_scale_box = NULL) {
  
  # If any names overlap between the aes_x and param_x lists, aes has higher priority
  # Beware of 'colour' and 'color' and what default aes() uses (if you for some reason supply colours in both aes and param)
  # Use new_col_scale and new_fill_scale if colour/fill are set inside aes() for both bee and box
  # and different colour scales are needed for the bees and the boxes
  # Should then supply col/fill_scale_bee (definitely) and col/fill_scale_box (can be defined after geom_beebox)
  
  # Make lists of inputs to the geom_x functions
  # Beeswarm plot parameters
  # Defaults are set here for user convenience, if set as a default argument,
  # all parameters would be overwritten when a list of just one parameter is supplied
  param_bee_default <- list(alpha = 1, size = 1.5, colour = "grey")
  # Overwrite defaults if any overlapping parameters are provided
  in_bee <- param_bee_default[setdiff(names(param_bee_default), names(param_bee))]
  # Add extra parameters
  in_bee <- append(in_bee, param_bee)
  # Parameters given in aes() are prioritised
  in_bee <- in_bee[setdiff(names(in_bee), names(aes_bee))]
  in_bee[["mapping"]] <- aes_bee
  
  # Boxplot parameters
  param_box_default <- list(alpha = 0, width = 0.1, linewidth = 1, colour = "black")
  in_box <- param_box_default[setdiff(names(param_box_default), names(param_box))]
  in_box <- append(in_box, param_box)
  in_box <- in_box[setdiff(names(in_box), names(aes_box))]
  in_box[["mapping"]] <- aes_box
  
  gg_out <- gg_in +
    list(do.call(geom_beeswarm, args = in_bee),
         if (!is.null(col_scale_bee)) col_scale_bee,
         if (!is.null(fill_scale_bee)) fill_scale_bee,
         if (new_col_scale) ggnewscale::new_scale_colour(),
         if (new_fill_scale) ggnewscale::new_scale_fill(),
         do.call(geom_boxplot, args = in_box),
         if (!is.null(col_scale_box)) col_scale_box,
         if (!is.null(fill_scale_box)) fill_scale_box)
  
  return(gg_out)
}


# Function for getting the loaded packages and package versions
pkg_version_tbl <- function() {
  pkg_names <- sort(names(sessionInfo()$otherPkgs))
  
  # Get package versions
  pkg_versions <- sapply(pkg_names, getNamespaceVersion)
  
  # Also get citations for each package
  pkg_citations <- sapply(pkg_names, function(x) {
    citation(x) %>% format(style = "text") %>%
      unlist() %>% paste(collapse = "")
  })
  
  tbl_out <- data.frame("Package" = pkg_names,
                        "Version" = pkg_versions,
                        "Citation" = pkg_citations)
  return(tbl_out)
}
