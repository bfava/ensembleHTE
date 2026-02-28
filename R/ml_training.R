
#' Encode factor/character columns as dummy variables
#' 
#' Used internally for algorithms that don't support factor features natively.
#' @param .X data.frame of covariates (may contain factor/character columns)
#' @return data.frame with factor/character columns replaced by dummy variables
#' @keywords internal
.encode_factors_if_needed <- function(.X) {
  factor_cols <- names(.X)[sapply(.X, function(col) is.factor(col) || is.character(col))]
  if (length(factor_cols) > 0) {
    # Convert character to factor first
    for (col in factor_cols) {
      if (is.character(.X[[col]])) {
        .X[[col]] <- as.factor(.X[[col]])
      }
    }
    .X <- fastDummies::dummy_cols(
      .X,
      select_columns = factor_cols,
      remove_selected_columns = TRUE,
      remove_first_dummy = TRUE
    )
  }
  .X
}

#' Fit a model using the appropriate framework
#' 
#' Dispatcher function that routes to mlr3 or grf based on algorithm
#' @param .Y Numeric vector of outcomes
#' @param .X data.frame of covariates
#' @param algorithm Character string specifying the ML algorithm
#' @param task_type Character: "regr" or "classif"
#' @param tune Logical, whether to tune hyperparameters
#' @param tune_params List of tuning parameters (see details in \code{ensemble_hte})
#' @param weights Optional numeric vector of observation weights
#' @param learner_params Optional named list of parameter-value pairs for this
#'   specific mlr3 learner (e.g., \code{list(num.trees = 500)} for ranger).
#'   Ignored for grf. Dispatched per-algorithm from \code{ensemble_hte}/\code{ensemble_pred}.
#' @keywords internal
fit_model <- function(.Y, .X, 
                      algorithm, 
                      task_type = c("regr", "classif"),
                      tune = FALSE, 
                      tune_params = list(
                        time = 30,
                        cv_folds = 3,
                        stagnation_iters = 250,
                        stagnation_threshold = 0.01,
                        measure = NULL
                      ),
                      weights = NULL,
                      learner_params = NULL) {
  
  if (algorithm == "grf") {
    .fit_model_grf(.Y, .X, task_type, weights)
  } else {
    .fit_model_mlr3(.Y, .X, algorithm, task_type, tune, tune_params, weights, learner_params)
  }
}

#' Fit a model using grf (generalized random forest)
#' @keywords internal
.fit_model_grf <- function(.Y, .X, 
                           task_type = c("regr", "classif"),
                           weights = NULL) {
  
  task_type <- match.arg(task_type)
  # grf requires numeric matrix — encode any factor/character columns
  .X <- .encode_factors_if_needed(.X)
  X_mat <- as.matrix(.X)
  
  if (task_type == "regr") {
    if (is.null(weights)) {
      model <- grf::regression_forest(X = X_mat, Y = .Y, num.threads = 1)
    } else {
      model <- grf::regression_forest(X = X_mat, Y = .Y, sample.weights = weights, num.threads = 1)
    }
  } else {
    # Classification: use probability_forest for binary outcomes
    Y_factor <- as.factor(.Y)
    if (is.null(weights)) {
      model <- grf::probability_forest(X = X_mat, Y = Y_factor, num.threads = 1)
    } else {
      model <- grf::probability_forest(X = X_mat, Y = Y_factor, sample.weights = weights, num.threads = 1)
    }
  }
  
  # Return a wrapper with consistent predict interface
  structure(
    list(
      model = model,
      task_type = task_type,
      levels = if (task_type == "classif") levels(as.factor(.Y)) else NULL
    ),
    class = "grf_wrapper"
  )
}

#' Predict method for grf_wrapper
#' @keywords internal
predict_newdata.grf_wrapper <- function(object, newdata) {
  newdata <- .encode_factors_if_needed(newdata)
  X_new <- as.matrix(newdata)
  
  if (object$task_type == "regr") {
    preds <- predict(object$model, newdata = X_new)$predictions
  } else {
    # For classification, return probability of positive class (second level)
    probs <- predict(object$model, newdata = X_new)$predictions
    preds <- probs[, 2]
  }
  
  list(response = preds)
}

#' Generic prediction helper for both mlr3 and grf models
#' @keywords internal
.predict_model <- function(model, newdata) {
  if (inherits(model, "grf_wrapper")) {
    predict_newdata.grf_wrapper(model, newdata)
  } else {
    # For mlr3 models, encode factors if the learner doesn't support them
    if (!("factor" %in% model$feature_types)) {
      newdata <- .encode_factors_if_needed(newdata)
    }
    pred <- model$predict_newdata(newdata)
    # For classification, extract probability of positive class (not factor response)
    if (inherits(pred, "PredictionClassif")) {
      # Get probability of second level (typically "1" for 0/1 outcomes)
      prob_matrix <- pred$prob
      pos_class <- colnames(prob_matrix)[2]  # Second column is positive class
      list(response = prob_matrix[, pos_class])
    } else {
      # For regression, just wrap in list for consistent interface
      list(response = pred$response)
    }
  }
}

#' Fit a model using mlr3
#' @param learner_params Optional named list of parameter-value pairs for this
#'   specific mlr3 learner (e.g., \code{list(num.trees = 500)} for ranger).
#'   Applied after algorithm-specific defaults and overrides them if there is a conflict.
#' @keywords internal
.fit_model_mlr3 <- function(.Y, .X, 
                            algorithm, 
                            task_type = c("regr", "classif"),
                            tune = FALSE, 
                            tune_params = list(
                              time = 30,
                              cv_folds = 3,
                              stagnation_iters = 250,
                              stagnation_threshold = 0.01,
                              measure = NULL
                            ),
                            weights = NULL,
                            learner_params = NULL) {
  
  # Validate task type
  task_type <- match.arg(task_type)
  
  # Handle NULL tune_params
  if (is.null(tune_params)) {
    tune_params <- list()
  }
  
  # Set default tuning parameters
  tune_params <- modifyList(
    list(
      time = 30,
      cv_folds = 3,
      stagnation_iters = 250,
      stagnation_threshold = 0.01,
      measure = NULL
    ),
    tune_params
  )

  # If "lm" is used with classification, treat as regression instead.
  # The user can specify "log_reg" explicitly if logistic regression is desired.
  if (task_type == "classif" && algorithm == "lm") {
    task_type <- "regr"
  }

  # Validate algorithm name before calling mlr3
  learner_id <- paste0(task_type, ".", algorithm)
  available_learners <- tryCatch(mlr_learners$keys(), error = function(e) character(0))
  if (length(available_learners) > 0 && !learner_id %in% available_learners) {
    stop("Algorithm '", algorithm, "' not recognized for task type '", task_type, "'. ",
         "Available algorithms include: lm, grf, ranger, glmnet, xgboost, nnet, kknn, svm. ",
         "For the full list of available learners, run: mlr3::mlr_learners$keys()")
  }

  # Encode factor columns only if the learner doesn't support them
  learner_tmp <- lrn(learner_id)
  if (!("factor" %in% learner_tmp$feature_types)) {
    .X <- .encode_factors_if_needed(.X)
  }

  # Prepare training data
  df_train <- data.frame(Y = .Y, .X)
  
  # Create task and learner based on type
  if (task_type == "regr") {
    task <- as_task_regr(df_train, target = "Y")
    learner <- lrn(learner_id)
    default_measure <- "regr.rsq"
  } else {
    task <- as_task_classif(df_train, target = "Y")
    learner <- lrn(learner_id)
    learner$predict_type <- "prob"  # Predict probabilities for binary outcomes
    default_measure <- "classif.auc"
  }
  
  # Add weights if provided
  if (!is.null(weights)) {
    task$cbind(data.frame(.weights = weights))
    task$set_col_roles(".weights", roles = "weights_learner")
  }
  
  # Set algorithm-specific parameters
  if (algorithm == 'nnet') {
    learner$param_set$values$trace <- FALSE
  }
  if (algorithm == 'gbm') {
    if (task_type == "regr") {
      learner$param_set$values$distribution <- 'gaussian'
    } else {
      learner$param_set$values$distribution <- 'bernoulli'
    }
  }
  if (algorithm == 'glmnet') {
    learner$param_set$values$s <- 0.01
  }
  if (algorithm == 'ranger') {
    learner$param_set$values$importance <- 'impurity'
  }
  
  # Apply user-specified learner parameters (override defaults if conflicting)
  if (!is.null(learner_params) && length(learner_params) > 0) {
    for (param_name in names(learner_params)) {
      tryCatch(
        learner$param_set$values[[param_name]] <- learner_params[[param_name]],
        error = function(e) {
          stop(paste0("Failed to set parameter '", param_name, "' on learner '",
                      learner_id, "'. ", conditionMessage(e),
                      "\nUse lrn('", learner_id, "')$param_set to see available parameters."))
        }
      )
    }
  }
  
  # Setup tuning if requested
  tuning_space_key <- paste0(task_type, '.', algorithm, '.default')
  if (tune && tuning_space_key %in% mlr3tuningspaces::mlr_tuning_spaces$keys()) {
    # Use user-provided mlr3tuning objects if available, otherwise use defaults
    tuner <- if (!is.null(tune_params$tuner)) {
      tune_params$tuner
    } else {
      mlr3tuning::tnr("random_search")
    }
    
    terminator <- if (!is.null(tune_params$terminator)) {
      tune_params$terminator
    } else {
      mlr3tuning::trm(
        "combo", 
        any = TRUE,
        list(
          mlr3tuning::trm("run_time", secs = tune_params$time %||% 30),
          mlr3tuning::trm(
            "stagnation", 
            iters = tune_params$stagnation_iters %||% 250, 
            threshold = tune_params$stagnation_threshold %||% 0.01
          )
        )
      )
    }
    
    resampling_method <- if (!is.null(tune_params$resampling)) {
      tune_params$resampling
    } else {
      rsmp("cv", folds = tune_params$cv_folds %||% 3)
    }
    
    # Use provided measure or default
    measure_name <- tune_params$measure %||% default_measure
    measure <- if (inherits(measure_name, "Measure")) measure_name else msr(measure_name)
    
    tuning_space <- if (!is.null(tune_params$search_space)) {
      tune_params$search_space
    } else {
      mlr3tuningspaces::lts(tuning_space_key)
    }
    
    model <- mlr3tuning::auto_tuner(
      tuner = tuner, 
      learner = learner, 
      resampling = resampling_method,
      measure = measure, 
      search_space = tuning_space, 
      terminator = terminator
    )
  } else {
    model <- learner
  }
  
  suppressWarnings(model$train(task))
  return(model)
}

#' Predict Individual Treatment Effects
#' 
#' @description
#' Internal function that predicts individual treatment effects (ITEs) using
#' the specified metalearner strategy and machine learning algorithm.
#' 
#' @param Y Numeric vector of outcomes
#' @param X data.frame of covariates
#' @param D Numeric vector of treatment indicators (0/1)
#' @param prop_score Numeric vector of propensity scores
#' @param W Numeric vector of inverse propensity weights
#' @param train_idx Logical vector indicating training observations
#' @param algorithm Character string specifying the ML algorithm
#' @param metalearner Character string specifying the metalearner strategy
#' @param r_learner Character string specifying the R-learner algorithm
#' @param task_type Character string: "regr" or "classif"
#' @param tune Logical, whether to tune hyperparameters
#' @param tune_params List of tuning parameters
#' @param test_idx Optional logical vector indicating which observations to
#'   predict for. If \code{NULL} (default), predicts for \code{!train_idx}.
#' @param learner_params Optional named list of parameter-value pairs for this
#'   specific mlr3 learner (e.g., \code{list(num.trees = 500)} for ranger).
#' 
#' @return List with predicted_y0, predicted_y1, and predicted_ite
#' 
#' @keywords internal
predict_ite <- function(Y, X, D, prop_score, W, train_idx, algorithm, 
                        metalearner = c("t", "s", "x", "r"), 
                        r_learner = "grf",
                        task_type = c("regr", "classif"), 
                        tune = FALSE,
                        tune_params = list(
                          time = 30,
                          cv_folds = 3,
                          stagnation_iters = 250,
                          stagnation_threshold = 0.01,
                          measure = NULL
                        ),
                        test_idx = NULL,
                        learner_params = NULL) {
  
  # Validate arguments
  task_type <- match.arg(task_type)
  metalearner <- match.arg(metalearner)
  
  # Convert character columns to factor (encoding is handled per-algorithm in fit_model)
  char_cols <- names(X)[sapply(X, is.character)]
  if (length(char_cols) > 0) {
    for (col in char_cols) {
      X[[col]] <- as.factor(X[[col]])
    }
  }
  
  # Split data
  X_train <- X[train_idx, , drop = FALSE]
  if (is.null(test_idx)) test_idx <- !train_idx
  X_test <- X[test_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  D_train <- D[train_idx]
  
  # Split train by treatment
  idx0 <- D_train == 0
  X_train0 <- X_train[idx0, , drop = FALSE]
  X_train1 <- X_train[!idx0, , drop = FALSE]
  Y_train0 <- Y_train[idx0]
  Y_train1 <- Y_train[!idx0]
  
  # Fit models + predict
  if (metalearner == "t") {
    model0 <- fit_model(Y_train0, X_train0, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    model1 <- fit_model(Y_train1, X_train1, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    predicted_y0 <- .predict_model(model0, X_test)$response
    predicted_y1 <- .predict_model(model1, X_test)$response
    predicted_ite <- predicted_y1 - predicted_y0

  } else if (metalearner == "s") {
    X_train_s <- cbind(X_train, D = D_train)
    model <- fit_model(Y_train, X_train_s, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    predicted_y0 <- .predict_model(model, cbind(X_test, D = 0))$response
    predicted_y1 <- .predict_model(model, cbind(X_test, D = 1))$response
    predicted_ite <- predicted_y1 - predicted_y0
    
  } else if (metalearner == "x") {
    # Step 1: Fit base learners
    model0 <- fit_model(Y_train0, X_train0, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    model1 <- fit_model(Y_train1, X_train1, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    
    # Step 2: Compute imputed treatment effects
    tau_train0 <- .predict_model(model1, X_train0)$response - Y_train0
    tau_train1 <- Y_train1 - .predict_model(model0, X_train1)$response
    
    # Step 3: Fit CATE models on imputed effects
    model_tau0 <- fit_model(tau_train0, X_train0, algorithm, "regr", tune, tune_params, learner_params = learner_params)
    model_tau1 <- fit_model(tau_train1, X_train1, algorithm, "regr", tune, tune_params, learner_params = learner_params)
    
    # Step 4: Combine predictions using propensity score
    prop_test <- prop_score[test_idx]
    predicted_y0 <- .predict_model(model0, X_test)$response
    predicted_y1 <- .predict_model(model1, X_test)$response
    predicted_ite <- prop_test * .predict_model(model_tau0, X_test)$response + 
      (1 - prop_test) * .predict_model(model_tau1, X_test)$response
    
  } else if (metalearner == "r") {
    # R-learner (Robinson transformation)
    
    # Step 1: Estimate E[Y|X] using the specified algorithm
    model_y <- fit_model(Y_train, X_train, algorithm, task_type, tune, tune_params, learner_params = learner_params)
    Y_hat_train <- .predict_model(model_y, X_train)$response
    
    # Step 2: Compute residuals (using actual propensity score, not estimated)
    prop_train <- prop_score[train_idx]
    resid_Y_train <- Y_train - Y_hat_train
    resid_D_train <- D_train - prop_train
    
    # Step 3: Estimate CATE
    if (r_learner == "grf") {
      # Use grf::causal_forest with pre-computed residuals
      # grf requires numeric matrix — encode any factor/character columns
      X_train_enc <- .encode_factors_if_needed(X_train)
      X_test_enc <- .encode_factors_if_needed(X_test)
      cf <- grf::causal_forest(
        X = as.matrix(X_train_enc),
        Y = Y_train,
        W = D_train,
        Y.hat = Y_hat_train,
        W.hat = prop_train,
        num.threads = 1
      )
      predicted_ite <- predict(cf, newdata = as.matrix(X_test_enc))$predictions
    } else {
      # Use mlr3 with weighted regression: resid_Y/resid_D ~ X, weights = resid_D^2
      pseudo_outcome <- resid_Y_train / resid_D_train
      weights_r <- resid_D_train^2
      
      model_tau <- fit_model(pseudo_outcome, X_train, r_learner, "regr", tune, tune_params, weights_r, learner_params = learner_params)
      predicted_ite <- .predict_model(model_tau, X_test)$response
    }
    
    # R-learner doesn't naturally produce y0/y1 predictions
    predicted_y0 <- NA
    predicted_y1 <- NA
  }
  
  list(
    predicted_y0 = predicted_y0, 
    predicted_y1 = predicted_y1, 
    predicted_ite = predicted_ite
  )
}
