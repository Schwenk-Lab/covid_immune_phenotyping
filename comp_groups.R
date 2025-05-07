comp_groups <- function(x, y, g = NA, 
                        ggtheme = theme_classic(),
                        cont_palette = "Paired", cont_legend = F,
                        discr_palette = "Set2",
                        tests = T, p_adj_method = "none", 
                        p_displ_thr = 0.05, p_digits = 2, p_stars = F, 
                        star_swap_list = list(c(0, 0.001, 0.01, 0.05, 1),
                                              c("***", "**", "*", "n.s")),
                        p_size = 6, write_p = T, sim_p = F, p_step = 0.07,
                        x_name = NULL, y_name = NULL, g_name = NULL, x_sizes = T) {
  ### --- Function for comparing groups --- ###
  ## Packages: ggplot2, ggsignif, dplyr (including the %>% pipe)
  # x, x-axis values or formula
  # y, y-axis values or data frame with data
  # g, additional grouping if desired. If x is a formula, g can be a character specifying the column name in input data frame y
  # ggtheme, theme to use for plot
  # cont_palette, RColorBrewer palette to use to fill boxes if y is continuous, or the scale_fill_x function call. NULL for no colour
  # cont_legend, logical indicating whether a legend should be drawn when the boxplot is filled in
  # discr_palette, RColorBrewer palette to use for groups in proportion bars if y is discrete, or the scale_fill_x function call
  # tests, if do statistical testing, currently not if multiple groupings
  # p_adj_method, test multiple adjustment
  # p_displ_thr, threshold for showing p values in plot
  # p_digits, number of digits to display for p-values
  # p_stars, logical, whether to display p-values as stars or not
  # star_swap_list, list containing as first element the thresholds for swapping p-values, and as second element the symbols to swap the numbers with
  # p_size, size of p-value annotation text
  # write_p, logical, if TRUE (default) "p = " is written before the actual p-value (unless showing stars instead) 
  # sim_p, logical, the argument 'simulate.p.value' of fisher.test()
  # p_step, numeric value for distance between p-value bars
  # x_name, y_name, g_name, names to use for axes and legends. If NULL, the names will be taken from the call
  
  ## Getting call for names of axes and title
  call_in <- match.call()
  
  # Determine whether extra grouping was specified or not. Extract x, y, and g arguments (split input if x is a formula)
  call_names <- call_in %>% as.list() %>% names()
  call_char <- as.character(call_in)
  
  if ("g" %in% call_names) {
    
    if (class(x) == "formula") {
      xsplit <- strsplit(call_char[2], "~")[[1]] %>% trimws()
      ax_names <- c(xsplit[2], xsplit[1], call_char[4])
    } else {
      ax_names <- call_char[2:4]
    }
    
    names(ax_names) <- c("x", "y", "g")
    if (!is.null(g_name)) {ax_names["g"] <- g_name}
    
  } else {
    
    if (class(x) == "formula") {
      xsplit <- strsplit(call_char[2], "~")[[1]] %>% trimws()
      ax_names <- c(xsplit[2], xsplit[1])
    } else {
      ax_names <- call_char[2:3]
    }
    
    names(ax_names) <- c("x", "y")
  }
  
  if (!is.null(x_name)) {ax_names["x"] <- x_name}
  if (!is.null(y_name)) {ax_names["y"] <- y_name}
  
  ## Make data frame for plotting
  if (class(x) == "formula") {
    
    # Allow for extra grouping variable specification via character (if formula input)
    if (length(g) == 1 & class(g) == "character") {
      ggdf <- data.frame(
        "x" = y[, as.character(x[3]), drop = T],
        "y" = y[, as.character(x[2]), drop = T],
        "g" = y[, g]
      )
    } else {
      ggdf <- data.frame(
        "x" = y[, as.character(x[3]), drop = T],
        "y" = y[, as.character(x[2]), drop = T],
        "g" = g
      )
    }
    
  } else {
    ggdf <- data.frame("x"=x, "y"=y, "g"=g)
  }
  
  ggdf <- ggdf %>% filter(!is.na(y))
  
  # Axis labels with number of values in each x group (if x_sizes == TRUE)
  x_unique <- as.character(sort(unique(ggdf$x)))
  x_labs <- if (x_sizes) {paste(x_unique, table(ggdf$x)[x_unique], sep = "\nn = ")} else {x_unique}
  names(x_labs) <- x_unique
  gg_xlabs <- scale_x_discrete(labels = x_labs)
  
  # x groups and combinations for tests
  # Make pairs after converting x values to characters as geom_signif does not work if the pairs are factors
  x_grp <- sort(unique(as.character(ggdf$x)))
  pairs <- combn(x_grp, 2, simplify = F)
  
  # g groups
  g_grp <- sort(unique(ggdf$g))
  
  # Colour scales
  # If a string is given, use it as the palette in `scale_fill_brewer`, otherwise use as is
  if (is.character(cont_palette)) {cont_palette = scale_fill_brewer(palette = cont_palette)}
  if (is.character(discr_palette)) {discr_palette = scale_fill_brewer(palette = discr_palette)}
  
  # Prepare a data frame for group combinations if extra group is specified, for statistical tests
  if (!all(is.na(ggdf$g))) {
    
    # Combinations of extra group values
    g_comb <- combn(g_grp, 2, simplify = F)
    
    # Combinations of x groups and extra groups
    x_grp_comb <- data.frame(matrix(nrow = (length(x_grp)*length(g_comb)),
                                    ncol = 4))
    colnames(x_grp_comb) <- c("x", "g1", "g2", "p.value")
    x_grp_comb$x <- rep(x_grp, each = length(g_comb))
    x_grp_comb$g1 <- rep(sapply(g_comb, function(i) i[1]), length(x_grp))
    x_grp_comb$g2 <- rep(sapply(g_comb, function(i) i[2]), length(x_grp))
    
  }
  
  ## -- Continuous y variable -- ##
  if (class(ggdf$y) %in% c("numeric", "double", "integer")) {
    
    # If grouping variable is given (or not all NA), group boxes by groups
    if (!all(is.na(ggdf$g))) {
      
      # Labels for facets
      facet_labs <- paste(ax_names["x"], x_grp, sep = ": ")
      names(facet_labs) <- x_grp
      
      # Turn g variable into factor for ensuring order is as desired (especially for p-value labels)
      ggdf$g <- factor(ggdf$g, levels = g_grp)
      
      # Make base plot, with or without colouring
      if (!is.null(cont_palette)) {
        plt <- ggplot(ggdf %>% filter(!is.na(g)), aes(x = g, y = y)) + 
          geom_boxplot(aes(fill = g), show.legend = cont_legend) +
          cont_palette +
          facet_wrap(~x, ncol = 5, labeller = labeller(x = facet_labs)) + 
          labs(x = ax_names["g"], y = ax_names["y"], fill = ax_names["g"]) + 
          ggtheme
      } else {
        plt <- ggplot(ggdf %>% filter(!is.na(g)), aes(x = g, y = y)) + 
          geom_boxplot() +
          facet_wrap(~x, ncol = 5, labeller = labeller(x = facet_labs)) + 
          labs(x = ax_names["g"], y = ax_names["y"]) + 
          ggtheme
      }
      
      ## -- Statistical tests -- ##
      if (tests & length(x_grp) > 1) {
        
        # Pairwise Wilcoxon test within x groups
        for (i in 1:nrow(x_grp_comb)) {
          x_grp_comb[i, "p.value"] <- 
            wilcox.test(y~g, 
                        ggdf %>% 
                          filter(x == x_grp_comb[i, "x"] &
                                   g %in% x_grp_comb[i, c("g1", "g2")]))$p.value
        }
        
        x_grp_comb$p.value <- p.adjust(x_grp_comb$p.value, method = p_adj_method)
        
        # Distance between bars as a percentage of total data range
        bardist <- (range(ggdf$y)[2] - range(ggdf$y)[1]) * 0.09
        
        # Max value in each facet to determine position of first bar
        facet_max <- ggdf %>% 
          group_by(x) %>% 
          summarise(ymax = max(y))
        
        # Make data frames for p-value labels and lines/bars
        
        # Convert group labels g1 and g2 to numbers to get their positions in plot x-axis
        p_lab_df <- x_grp_comb
        p_lab_df$g1 <- as.numeric(factor(p_lab_df$g1, levels = g_grp))
        p_lab_df$g2 <- as.numeric(factor(p_lab_df$g2, levels = g_grp))
        
        # Compute coordinate for p-value label (right in-between)
        p_lab_df$g <- (p_lab_df$g1 + p_lab_df$g2)/2
        
        # Add p-value star representations
        p_lab_df$p.star <- symnum(p_lab_df$p.value, 
                                  cutpoints = star_swap_list[[1]],
                                  symbols = star_swap_list[[2]],
                                  corr = F, legend = F) %>% as.character()
        
        # Remove p-values that will not be displayed
        p_lab_df <- p_lab_df %>% filter(p.value < p_displ_thr)
        
        # Compute y coordinates of p-values
        p_lab_df$y <- sapply(unique(p_lab_df$x), function(x) {
          y_bar1 <- as.numeric(facet_max[facet_max$x == x, "ymax"]) + bardist
          y_bar_max <- y_bar1 + bardist * (length(which(p_lab_df$x == x)) - 1)
          seq(y_bar1, y_bar_max, by = bardist)
        }) %>% unlist() %>% as.numeric()
        
        if (nrow(p_lab_df) > 0) {
          # Coordinates for lines
          p_line_df <- apply(p_lab_df, 1, function(rowx) {
            data.frame("x" = rep(rowx["x"], 2), 
                       "g" = as.numeric(rowx[c("g1", "g2")]),
                       "y" = rep(as.numeric(rowx["y"]) - (bardist / 3), 2),
                       "grp" = paste0(rowx["x"], rowx["g1"], rowx["g2"]))
          }) %>% bind_rows() 
          
          p_lab_df$p.value <- signif(p_lab_df$p.value, p_digits)
          
          plt <- plt + 
            geom_text(data = p_lab_df, aes_string(label = ifelse(p_stars, "p.star", "p.value"))) + 
            geom_line(data = p_line_df, aes(group = grp))
        }
        
        return(list("p.value"=x_grp_comb,
                    "plt"=plt + ggtheme))
      }
      
    } else {
      
      # Make base plot
      if (!is.null(cont_palette)) {
        plt <- ggplot(ggdf, aes(x = x, y = y)) + 
          geom_boxplot(aes(fill = x), show.legend = cont_legend) +
          cont_palette +
          labs(x = ax_names["x"], y = ax_names["y"], fill = ax_names["x"]) +
          gg_xlabs + ggtheme
      } else {
        plt <- ggplot(ggdf, aes(x = x, y = y)) + 
          geom_boxplot() +
          labs(x = ax_names["x"], y = ax_names["y"]) + 
          gg_xlabs + ggtheme
      }
      
      ## -- Statistical tests -- ##
      if (tests & length(x_grp) > 1) {
        
        # Pairwise Wilcoxon test
        pvals <- vector("numeric", length(pairs))
        for (i in 1:length(pairs)) {
          pvals[i] <- wilcox.test(
            x = ggdf %>% filter(x == pairs[[i]][1]) %>% pull(y), 
            y = ggdf %>% filter(x == pairs[[i]][2]) %>% pull(y)
          )$p.value
        }
        
        # Multiple testing adjustment
        pvals <- p.adjust(pvals, method = p_adj_method)
        
        # P-values for printing
        pvals_out <- matrix(nrow = length(x_grp), 
                            ncol = length(x_grp),
                            dimnames = list(x_grp, x_grp))
        for (i in seq(1, length(pairs))) {
          pvals_out[pairs[[c(i, 2)]], pairs[[c(i, 1)]]] <- 
            pvals[i]
        }
        pvals_out <- pvals_out[-1, -ncol(pvals_out)]
        
        # Star representation of p-values
        pstars <- symnum(pvals,
                         cutpoints = star_swap_list[[1]],
                         symbols = star_swap_list[[2]],
                         corr = F, legend = F) %>% as.character()
        
        # Plot those below given threshold
        pairs <- pairs[which(pvals < p_displ_thr)]
        pstars <- pstars[which(pvals < p_displ_thr)]
        pvals <- pvals[which(pvals < p_displ_thr)] %>% signif(p_digits)
        
        # Add p-values to plot if any are left
        if (length(pvals) > 0) {
          if (write_p) {
            pvals <- paste0("italic(p) == ", pvals)
          }
          
          plt <- plt + 
            geom_signif(comparisons = pairs, 
                        step_increase = p_step, 
                        annotations = ifelse(rep(p_stars, length(pairs)), pstars, as.character(pvals)),
                        parse = !p_stars, textsize = p_size)
        }
        
        return(list("p.value"=pvals_out, 
                    "plt"=plt + ggtheme))
        
      }
    }
  }
  
  ## -- Categorical y variable -- ##
  if (class(ggdf$y) %in% c("character", "factor")) {
    
    # Multiple groupings if specified
    if (!all(is.na(ggdf$g))) {
      
      # Grouped data frame that can be used for tests?
      ggdf_grouped <- ggdf %>%
        filter(!is.na(y) & !is.na(g)) %>%
        group_by(x, g) %>%
        count(y) %>%
        mutate(freq = n / sum(n))
      
      # Facet labels
      facet_labs <- paste(ax_names["x"], x_grp, sep = ": ")
      names(facet_labs) <- x_grp
      
      plt <- ggplot(ggdf_grouped, aes(x = g, y = freq, fill = as.factor(y))) + 
        geom_col(position = "fill") + 
        facet_wrap(~x, ncol = 5, labeller = labeller(x = facet_labs)) + 
        labs(x = ax_names["g"], y = "Proportion", fill = ax_names["y"]) + 
        ggtheme
      
      ## -- Statistical tests -- ##
      if (tests & length(x_grp) > 1) {
        
        # Pairwise Fisher exact test
        for (i in 1:nrow(x_grp_comb)) {
          
          # Make contigency table for relevant x and extra group combinations
          cont_tbl <- table(
            ggdf %>% filter(x == x_grp_comb[i, "x"] & g %in% x_grp_comb[i, c("g1", "g2")]) %>% pull(y), 
            ggdf %>% filter(x == x_grp_comb[i, "x"] & g %in% x_grp_comb[i, c("g1", "g2")]) %>% pull(g) 
          )
          
          x_grp_comb[i, "p.value"] <- 
            fisher.test(cont_tbl, simulate.p.value = sim_p)$p.value
        }
        
        x_grp_comb$p.value <- p.adjust(x_grp_comb$p.value, method = p_adj_method)
        
        # Vertical distance between bars, plot goes between 0 and 1
        bardist <- 0.07
        
        # Data frame for p-value labels
        p_lab_df <- x_grp_comb
        p_lab_df$g1 <- as.numeric(factor(p_lab_df$g1, levels = g_grp))
        p_lab_df$g2 <- as.numeric(factor(p_lab_df$g2, levels = g_grp))
        
        # Compute coordinate for p-value label (right in-between)
        p_lab_df$g <- (p_lab_df$g1 + p_lab_df$g2)/2
        
        # Star representations of p-values
        p_lab_df$p.star <- symnum(p_lab_df$p.value, 
                                  cutpoints = star_swap_list[[1]],
                                  symbols = star_swap_list[[2]],
                                  corr = F, legend = F) %>% as.character()
        
        # Remove p-values that will not be displayed
        p_lab_df <- p_lab_df %>% filter(p.value < p_displ_thr)
        
        # Compute y coordinates of p-values
        p_lab_df$y_pos <- sapply(unique(p_lab_df$x), function(x) {
          y_bar1 <- 1 + bardist
          y_bar_max <- y_bar1 + bardist * (length(which(p_lab_df$x == x)) - 1)
          seq(y_bar1, y_bar_max, by = bardist)
        }) %>% unlist() %>% as.numeric()
        
        if (nrow(p_lab_df) > 0) {
          # Coordinates for lines
          p_line_df <- apply(p_lab_df, 1, function(rowx) {
            data.frame("x" = rep(rowx["x"], 2), 
                       "g" = as.numeric(rowx[c("g1", "g2")]),
                       "y_pos" = rep(as.numeric(rowx["y_pos"]) - (bardist / 3), 2),
                       "grp" = paste0(rowx["x"], rowx["g1"], rowx["g2"]))
          }) %>% data.table::rbindlist()
          
          p_lab_df$p.value <- signif(p_lab_df$p.value, p_digits)
          plt <- plt +
            geom_text(data = p_lab_df, aes_string(x = "g", y = "y_pos", fill = NULL, label = ifelse(p_stars, "p.star", "p.value"))) +
            geom_line(data = p_line_df, aes(x = g, y = y_pos, fill = NULL, group = grp))
        }
        
        colnames(x_grp_comb) <- c(ax_names["x"], 
                                  paste0(ax_names["g"], "_1"), 
                                  paste0(ax_names["g"], "_2"), 
                                  "p.value")
        return(list("p.value" = x_grp_comb,
                    "plt" = plt + discr_palette +
                      ggtheme + scale_y_continuous(breaks = seq(0, 1, by = 0.25))))
      }
      
    } else {
      
      # Compute proportions manually to be compatible with geom_signif later
      plt <- ggplot(ggdf %>% 
                      filter(!is.na(y)) %>%
                      group_by(x) %>%
                      count(y) %>%
                      mutate(freq = n / sum(n)), 
                    aes(x = x, y = freq, fill = y)) + 
        geom_col(position = "fill") + 
        labs(x = ax_names["x"], y = "Proportion", fill = ax_names["y"]) + 
        gg_xlabs + ggtheme
      
      ## -- Statistical tests -- ##
      if (tests & length(x_grp) > 1) {
        
        # Pairwise Fisher exact test
        pvals <- vector("numeric", length(pairs))
        
        # Contingency table
        cont_tbl <- table(ggdf$x, ggdf$y)
        
        for (i in 1:length(pairs)) {
          tbl <- cont_tbl[pairs[[i]], ]
          pvals[i] <- fisher.test(tbl, simulate.p.value = sim_p)$p.value
        }
        
        # Multiple testing adjustment
        pvals <- p.adjust(pvals, method = p_adj_method)
        
        # P-values for printing
        pvals_out <- matrix(nrow = length(x_grp), 
                            ncol = length(x_grp),
                            dimnames = list(x_grp, x_grp))
        for (i in seq(1, length(pairs))) {
          pvals_out[pairs[[c(i, 2)]], pairs[[c(i, 1)]]] <- 
            pvals[i]
        }
        pvals_out <- pvals_out[-1, -ncol(pvals_out)]
        
        # Star representations of p-values
        pstars <- symnum(pvals,
                         cutpoints = star_swap_list[[1]],
                         symbols = star_swap_list[[2]],
                         corr = F, legend = F) %>% as.character()
        
        # Only display if above given threshold
        pairs <- pairs[which(pvals < p_displ_thr)]
        pstars <- pstars[which(pvals < p_displ_thr)]
        pvals <- pvals[which(pvals < p_displ_thr)] %>% signif(p_digits)
        
        # Add p-values to plot if any left
        if (length(pvals) > 0) {
          if (write_p) {
            pvals <- paste0("italic(p) == ", pvals)
          }
          
          plt <- plt + 
            geom_signif(comparisons = pairs, 
                        annotations = ifelse(rep(p_stars, length(pairs)), pstars, as.character(pvals)),
                        parse = !p_stars, textsize = p_size,
                        y_position = seq(1.03, 1.03 + length(pvals)*p_step, by = p_step))
        }
        
        return(list("p.value" = pvals_out, 
                    "plt" = plt + discr_palette +
                      ggtheme + scale_y_continuous(breaks = seq(0, 1, by = 0.25))))
        
      }
    }
  }
  
  return(plt + ggtheme)
}




