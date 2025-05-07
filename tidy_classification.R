#' Wrapper function for ML classification using the tidymodels package (lasso or xgboost)
#'
#' @param x Data frame or tibble with predictors (columns of any classes) and response (character or factor column).
#' @param response The unquoted name of the response column.
#' @param train_test_split Number between 0 and 1 specifying the fraction of samples to be used for training and tuning. Alternatively a vector of integers or doubles indicating the indices of the training set samples.
#' @param strat_split Whether the split should be stratified by the outcome or not.
#' @param engine Name of engine to use (default is `glmnet`).
#' @param method Character of length one specifying method to use (Lasso or XGBoost).
#' @param tune_metrics A call to the `yardstick::metric_set()` function with the `yardstick` metric functions to use for evaluating tuning results
#' @param tune_select Either a character of length one with the metric to use for picking the optimal parameter from tuning, or a single numeric value specifying the index of the model to use (for picking a specific penalty value based on a previous run with the same seed).
#' @param preproc_steps A sequence of preprocessing step functions of the `recipes` package, piped together by `magrittr` pipes. The default behaviour is to filter out zero-variance predictors, create dummy variables of nominal predictors, and scale and center numeric predictors. Preprocessing steps can be skipped by giving a pipeline with an anonymous function returning its input (i.e. `. %>% (\(x) x)()`) (move this last part to examples?). Also allows for updating roles of variables by including `update_role()` into the pipe.
#' @param tune_v Number of folds to use for tuning the penalty.
#' @param tune_grid_size Number of values to test for tuning.
#' @param rseed Random seed to use, set to a non-number to not set any seed.
#'
#' @return A list containing model metrics, a ROC curve, variable importance, tuning results, plot for tuning, the preprocessing steps, and the fit object.
#'
#' @examples tidy_classification(iris, Species)
tidy_classification <- function(x, response, train_test_split = 3/4, strat_split = T,
                                engine = "glmnet", method = c("lasso", "xgboost"),
                                tune_metrics = metric_set(roc_auc), tune_select = "roc_auc",
                                preproc_steps = . %>%
                                  step_zv(all_predictors()) %>%
                                  step_dummy(all_nominal_predictors()) %>%
                                  step_normalize(all_numeric_predictors()),
                                tune_v = 5, tune_grid_size = 30, rseed = 123) {
  
  # Name of response column for use where a string is needed
  response_name <- rlang::as_name(enquo(response))
  
  # Check if response is binomial or multinomial
  if (length(unique(pull(x, {{response}}))) == 2) {
    class_nom <- "binomial"
  } else if (length(unique(pull(x, {{response}}))) > 2) {
    class_nom <- "multinomial"
  } else {
    stop("Response needs to have at least two classes!")
  }
  
  # Set seed
  if (is.numeric(rseed)) {
    set.seed(rseed)
  }
  
  # Define model
  if (method[1] == "lasso" & class_nom == "binomial") {
    model <- logistic_reg(mode = "classification", engine = engine, penalty = tune(), mixture = 1)
    tune_params <- "penalty"
    
  } else if (method[1] == "lasso" & class_nom == "multinomial") {
    model <- multinom_reg(mode = "classification", engine = engine, penalty = tune(), mixture = 1)
    tune_params <- "penalty"
    
  } else if (method[1] == "xgboost") {
    model <- boost_tree(mode = "classification", engine = "xgboost", trees = 1001,
                        min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune())
    tune_params <- c("min_n", "tree_depth", "learn_rate", "loss_reduction")
    
  } else {
    stop("`method` must be 'lasso' or 'xgboost'.")
  }
  
  # Make grid for tuning
  if (method[1] == "xgboost") {
    tune_grid <- grid_max_entropy(list(min_n(), tree_depth(), learn_rate(), loss_reduction()), size = tune_grid_size)
  } else {
    tune_grid <- grid_regular(penalty(), levels = tune_grid_size)
  }
  
  # Make training data for training and tuning
  if (length(train_test_split) == 1) {
    if (strat_split) {
      splits <- initial_split(x, prop = train_test_split, strata = {{response}})
    } else {
      splits <- initial_split(x, prop = train_test_split)
    }
    
    dat_train <- training(splits)
    
  } else if (length(train_test_split) > 1) {
    # Make split object with arbitrary proportion and overwrite indices in resulting object
    splits <- initial_split(x, prop = 0.5)
    splits$in_id <- train_test_split
    dat_train <- training(splits)
    
  } else {
    stop("Not a valid split!")
  }
  
  # Make recipe
  rcp <- recipe(as.formula(paste0(response_name, " ~ .")), data = dat_train) %>%
    preproc_steps
  
  # Make workflow
  wf <- workflow() %>%
    add_model(model) %>%
    add_recipe(rcp)
  
  # Tune model
  tune_set <- vfold_cv(dat_train, v = tune_v)
  tune_results <- wf %>%
    tune_grid(resamples = tune_set, grid = tune_grid,
              control = control_grid(save_pred = T),
              # Set metrics from input list
              metrics = tune_metrics)
  
  # Get model with best (or selected) parameter
  tune_best <- if (is.character(tune_select)) {
    # If a metric is given, pick the best model for that metric
    # Selects the first model that has the highest value of the selected metric
    # select_best(tune_results, tune_select)
    collect_metrics(tune_results) %>%
      filter(.metric == tune_select) %>%
      slice_max(mean, n = 1, with_ties = F)
    
  } else if (is.numeric(tune_select)) {
    # If a number is given, pick the model with that index
    # Add a leading zero if the index is smaller than 10
    tune_results %>%
      collect_metrics() %>%
      filter(str_detect(.config, paste0(ifelse(nchar(tune_select) == 1,
                                               paste0("0", tune_select),
                                               tune_select),
                                        "$"))) %>%
      dplyr::slice(1)
    
  }
  
  # Plot AUC against tuned parameters
  # Different depending on method used
  # Make common data
  tune_plot_dat <- tune_results %>%
    collect_metrics() %>%
    select(all_of(tune_params), .metric, .config, mean)
  
  used_metrics <- unique(c(tune_select, tune_plot_dat$.metric))
  tune_plot <- tune_plot_dat %>%
    pivot_longer(cols = -c(.config, .metric)) %>%
    mutate(name = case_when(name == "mean" ~ .metric, T ~ name),
           name = factor(name, levels = c(used_metrics, tune_params)),
           model = as.numeric(str_extract(.config, "\\d+$"))) %>%
    select(-.metric) %>%
    ggplot(aes(x = model, y = value)) +
    geom_point() + geom_line(aes(group = name)) +
    geom_vline(xintercept = as.numeric(str_extract(tune_best$.config, "\\d+$")),
               linetype = "dashed", colour = "grey") +
    facet_wrap(~ name, ncol = 1, scales = "free") +
    theme_classic()
  
  # Make the final workflow
  wf_final <- wf %>% finalize_workflow(dplyr::slice(tune_best, 1))
  
  # Fit and evaluate the model
  final_fit <- wf_final %>% last_fit(splits)
  
  # Make ROC curve
  if (class_nom == "binomial") {
    pred_class <- if (is.character(pull(x, {{response}}))) {
      # If character, take first value when sorted as predicted class, to be consistent with the assumption of the roc_curve() function
      paste0(".pred_", sort(unique(pull(x, {{response}})))[1])
    } else {
      # If factor, take first factor level
      paste0(".pred_", levels(pull(x, {{response}}))[1])
    }
    
    roc_curve <- final_fit %>%
      collect_predictions() %>%
      roc_curve({{response}}, !!pred_class) %>%
      ggplot(aes(x = 1 - specificity, y = sensitivity)) +
      geom_path(linewidth = 1) + geom_abline(slope = 1, colour = "grey", linetype = "dashed") +
      labs(colour = response_name, y = "Sensitivity") +
      theme_classic(16) +
      theme(aspect.ratio = 1)
    
  } else if (class_nom == "multinomial") {
    
    # Get per-class AUC
    auc_per_class <- final_fit %>%
      # Get probabilities
      collect_predictions() %>%
      select({{response}}, matches("\\.pred_"), -.pred_class) %>%
      pivot_longer(cols = -{{response}}) %>%
      # Make names match actual labels
      mutate(name = str_remove(name, "\\.pred_")) %>%
      # Calculate AUC per class
      mutate(auc = roc_auc_vec(
        truth = factor(get(response_name) == name, levels = c(F, T)),
        estimate = value, event_level = "second"
      ), .by = name) %>%
      distinct(name, auc)
    
    # Make ROC curves, add computed per-class AUCs
    roc_curve <- final_fit %>%
      collect_predictions() %>%
      roc_curve({{response}}, matches("\\.pred"), -.pred_class) %>%
      left_join(auc_per_class, by = c(".level" = "name")) %>%
      mutate(.level = paste0(.level, ", AUC: ", signif(auc, 3))) %>%
      ggplot(aes(x = 1 - specificity, y = sensitivity, colour = .level)) +
      geom_path(linewidth = 1) + geom_abline(slope = 1, colour = "grey", linetype = "dashed") +
      scale_colour_brewer(palette = "Pastel1") +
      labs(colour = response_name, y = "Sensitivity") +
      theme_classic(16) +
      theme(aspect.ratio = 1)
    
  }
  
  # Get variable importance for all variables
  varimp <- {if (method[1] == "lasso") {
    # Lasso models
    final_fit %>%
      extract_fit_parsnip() %>%
      tidy() %>%
      filter(term != "(Intercept)") %>%
      select(-all_of(tune_params)) %>%
      rename(variable = term, importance = estimate) %>%
      # Convert to absolute values
      mutate(sign = case_when(importance < 0 ~ "Negative",
                              importance > 0 ~ "Positive",
                              T ~ "0"),
             importance = abs(importance))
    
  } else if (method[1] == "xgboost") {
    # Boosted tree model
    final_fit %>%
      extract_fit_engine() %>%
      xgb.importance(model = .) %>%
      rename(variable = Feature, importance = Gain, cover = Cover, frequency = Frequency) %>%
      mutate(sign = 0)
    
  }} %>%
    # Add a dummy class variable if binomial outcome to work with scaling code, or if importances are not obtained per class (xgboost)
    mutate(class = if (class_nom == "binomial" | method[1] == "xgboost") {"dummy"} else {class}) %>%
    # Add scaled importance, ranging between 0 and 1 in each class
    # To make sure 0 is 0 in the new scale, add a row with importance 0 and remove it after scaling
    group_by(class) %>% nest() %>%
    mutate(importance_scaled = map(data, \(x) {
      x %>% add_row(importance = 0) %>%
        mutate(imp_sc = 1 * ((importance - min(importance)) / (max(importance) - min(importance)))) %>%
        filter(!is.na(variable)) %>%
        pull(imp_sc)
    })) %>%
    unnest(c(data, importance_scaled)) %>%
    ungroup()
  
  # If there were only 0 importances to begin with the scaling will return NaN, turn those back into 0  
  varimp$importance_scaled[is.nan(varimp$importance_scaled)] <- 0
  
  # Make variable importance plot with top 10
  varimp_plot <- varimp %>%
    group_by(class) %>% nest() %>%
    mutate(plt = map(data, \(x) {
      x %>% filter(importance_scaled > 0) %>%
        slice_max(importance_scaled, n = 10) %>%
        mutate(variable = factor(variable, levels = variable),
               sign = factor(sign, levels = c("Positive", "Negative"))) %>%
        # Order by increasing importance to start at top
        ggplot(aes(x = importance_scaled, y = reorder(variable, importance_scaled), colour = sign)) +
        geom_segment(aes(yend = variable, x = 0, xend = importance_scaled),
                     linewidth = 1, show.legend = T) +
        geom_point(size = 3, show.legend = T) +
        scale_colour_manual(values = c("Positive" = "sienna2", "Negative" = "skyblue2"),
                            drop = F) +
        labs(title = ifelse(class_nom == "binomial", "", paste0(response_name, ": ", class)),
             y = "Variable", x = "Scaled importance", colour = "Sign") +
        xlim(0, 1) +
        theme_classic(16) +
        theme(plot.title = element_text(hjust = 0.5))
    })) %>% pull(plt)
  
  # Make a confusion matrix plot
  confmat <- final_fit %>%
    collect_predictions() %>%
    conf_mat({{response}}, .pred_class) %>%
    pluck("table") %>% as_tibble()
  
  confmat_plot <- ggplot(confmat, aes(x = Truth, y = Prediction, fill = n)) +
    geom_tile() + geom_text(aes(label = n)) +
    scale_fill_gradient(low = "white", high = "sienna2") +
    coord_fixed() +
    theme_classic()
  
  # Remove dummy class variable if binomial outcome
  if (class_nom == "binomial" | method[1] == "xgboost") {varimp <- select(varimp, -class)}
  
  return(list(
    "metrics" = final_fit %>% collect_metrics(),
    "roc_curve" = roc_curve,
    "varimp" = varimp,
    "confmat" = confmat,
    "vip" = varimp_plot,
    "conf_plot" = confmat_plot,
    "fit" = final_fit,
    "tune_selected" = tune_best,
    "tune_plot" = tune_plot,
    "preprocessing_steps" = preproc_steps
  ))
}
