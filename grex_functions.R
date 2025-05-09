# grex_functions.R
# functions for GReX step of pipeline
# Nolan Cole, UW Biostatistics 26 march 2025

# grex_models(ncrna = "ENSG00000241599", Y = expression_data, X = snp_data, dosages = snp_dosage, sample_ids = sample_ids)
grex_models <- function(ncrna, Y, X, dosages, sample_ids) {
  
  ## Fit predictive models
  # SuSiE
  susie <- univariate_susie(X,
                            Y,
                            Omega = 1,
                            scale = FALSE,
                            alpha = 0.5,
                            nfolds = 10,
                            verbose = FALSE,
                            par = FALSE,
                            n.cores = NULL,
                            tx_names = ncrna,
                            seed = 2025)[[1]]
  
  if (ncol(X) > 1) {
    # Fit Lasso (alpha = 1)
    lasso <- cv.glmnet(X, Y, alpha = 1, nfolds = 10, keep = TRUE)
    # Adjusted R^2
    lambda_idx_lasso <- which(lasso$lambda == lasso$lambda.min)
    y_hat_cv_lasso <- lasso$fit.preval[, lambda_idx_lasso]
    r2_cv_lasso <- 1 - sum((Y - y_hat_cv_lasso)^2) / sum((Y - mean(Y))^2)
    p_lasso <- lasso$nzero[lambda_idx_lasso]
    n <- length(Y)
    lasso_adj_R2 <- 1 - ((1 - r2_cv_lasso) * (n - 1) / (n - p_lasso - 1))
    names(lasso_adj_R2) <- NULL
    # check if all zero coefficients
    lasso_all_zero <- all(as.vector(coef(lasso, s = "lambda.min"))[-1] == 0)
    
    # Fit Elastic Net (alpha = 0.5)
    enet <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10, keep = TRUE)
    # Adjusted R^2
    lambda_idx_enet <- which(enet$lambda == enet$lambda.min)
    y_hat_cv_enet <- enet$fit.preval[, lambda_idx_enet]
    r2_cv_enet <- 1 - sum((Y - y_hat_cv_enet)^2) / sum((Y - mean(Y))^2)
    p_enet <- enet$nzero[lambda_idx_enet]
    enet_adj_R2 <- 1 - ((1 - r2_cv_enet) * (n - 1) / (n - p_enet - 1))
    names(enet_adj_R2) <- NULL
    
    # HAL
    hal <- hal9001::fit_hal(X, Y,
                            smoothness_orders = 1, 
                            max_degree = 1,
                            num_knots = 5,
                            fit_control = list(nfolds = 10),
                            return_lasso = TRUE)
    hal_R2_cv <- cross_fitted_r2(Y, X, hal$lambda_star)
    
    
  } else if (ncol(X) == 1) {
    
    lasso_adj_R2 = enet_adj_R2 = hal_R2_cv = 0
    lasso_all_zero = enet_all_zero <- TRUE
    
  }
  
  
  r2vec <- c("lasso" = lasso_adj_R2,
             "enet" = enet_adj_R2,
             "susie" = susie$R2,
             "hal" = hal_R2_cv)
  mmm <- r2vec[which.max(r2vec)] |> names()
  
  if (max(r2vec) > 0.01) {
    # just use susie if 0 weights or HAL
    if(mmm == "susie" | mmm == "hal" | (mmm == "lasso" & lasso_all_zero)) {
      best_model <- susie
      row <- c(ncrna, mmm, scale(c(susie$Pred)))
    } else if (mmm == "enet") {
      # Best Lambda (minimum MSE)
      best_lambda <- enet$lambda.min
      # Final lasso model with best lambda
      best_model <- glmnet(X, Y,
                           alpha = 0.5,
                           lambda = best_lambda)
      enet_pred <- predict(best_model, X)
      row <- c(ncrna, mmm, enet_pred)
    } else if (mmm == "lasso") {
      # Best Lambda (minimum MSE)
      best_lambda <- lasso$lambda.min
      # Final lasso model with best lambda
      best_model <- glmnet(X, Y,
                           alpha = 1,
                           lambda = best_lambda)
      lasso_pred <- predict(best_model, X)
      row <- c(ncrna, mmm, lasso_pred)
    }
    row <- as.matrix(row)
    row <- as.data.frame(t(row))
    colnames(row) <- c("gene", "model", sample_ids)
    
    return(list("row" = row, "best_model" = best_model))
  } else {
    return(NULL)
  }
}

cross_fitted_r2 <- function(Y, X, min_lambda = NULL) {
  
  # Compute adjusted R^2 using cross-fitting
  
  K <- 10
  n <- nrow(X)
  folds <- sample(rep(1:K, length.out = n))
  Y_hat_cv <- rep(NA, n)
  
  for (k in 1:K) {
    train_idx <- which(folds != k)
    val_idx <- which(folds == k)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    X_val <- X[val_idx, , drop = FALSE]
    hal_fit_k <- fit_hal(X = X_train, Y = y_train,
                         lambda = min_lambda,
                         smoothness_orders = 1, 
                         max_degree = 1,
                         num_knots = 5,
                         fit_control = list(cv_select = FALSE))
    
    Y_hat_cv[val_idx] <- predict(hal_fit_k, 
                                 new_data = X_val, 
                                 s = min_lambda)
  }
  
  # Compute cross-fitted R^2
  rss <- sum((Y - Y_hat_cv)^2)
  tss <- sum((Y - mean(Y))^2)
  R2_cv <- 1 - rss / tss
  
  return(R2_cv)
  
}
