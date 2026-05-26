# 01_ps_estimators.R
# Propensity score estimation methods for sequential trial emulation
# Implements 8 estimators: MLE, BC-MLE, FL, FLIC, FLAC, LF, CP, LASSO

library(dplyr)
library(brglm2)   # Firth, FLIC, FLAC
library(logistf)  # log-F(1,1)
library(arm)      # bayesglm (Cauchy)
library(glmnet)   # LASSO
library(boot)     # Bootstrap

#' Estimate propensity scores and weights: Maximum Likelihood (MLE)
#' 
#' Standard logistic regression via maximum likelihood.
#' Reference estimator; diverges under separation.
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model (e.g., A_base ~ Lc_base + Lb_base + V)
#' @return List with ps (fitted probabilities), weights (SIPTCW numerator component),
#'   converged (logical), method name
#' @export
ps_mle <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  fit <- tryCatch(
    glm(formula, data = data, family = binomial()),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  if (is.null(fit) || !fit$converged) {
    return(list(
      ps = rep(NA, nrow(data)),
      weights = rep(NA, nrow(data)),
      converged = FALSE,
      method = "MLE"
    ))
  }
  
  ps <- predict(fit, type = "response")
  # Bound away from 0/1 to avoid infinite weights
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  # Stabilized treatment weight numerator: marginal prevalence
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  
  # sw^A = P(A=a) / P(A=a | L, V)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = ps,
    weights = weights,
    converged = TRUE,
    method = "MLE"
  ))
}


#' Estimate propensity scores: Bootstrap Bias-Corrected MLE (BC-MLE)
#' 
#' Implements Sanim et al. (2026) bootstrap bias correction:
#' beta_bc = 2 * beta_mle - mean(beta_bootstrap)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @param B Number of bootstrap resamples (default 200)
#' @return List with ps, weights, converged, method
#' @export
ps_bc_mle <- function(data, formula = A_base ~ Lc_base + Lb_base + V, B = 200) {
  
  # Fit MLE
  fit_mle <- tryCatch(
    glm(formula, data = data, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit_mle) || !fit_mle$converged) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "BC-MLE"))
  }
  
  beta_mle <- coef(fit_mle)
  
  # Bootstrap at patient level (to preserve within-patient correlation)
  ids <- unique(data$id)
  n_ids <- length(ids)
  
  beta_boot <- matrix(NA, nrow = B, ncol = length(beta_mle))
  
  for (b in 1:B) {
    # Resample IDs
    ids_b <- sample(ids, size = n_ids, replace = TRUE)
    data_b <- data %>%
      filter(id %in% ids_b) %>%
      group_by(id) %>%
      mutate(id_new = cur_group_id()) %>%
      ungroup()
    
    fit_b <- tryCatch(
      glm(formula, data = data_b, family = binomial()),
      error = function(e) NULL
    )
    
    if (!is.null(fit_b) && fit_b$converged) {
      beta_boot[b, ] <- coef(fit_b)
    }
  }
  
  # Bias-corrected coefficients
  beta_b_mean <- colMeans(beta_boot, na.rm = TRUE)
  beta_bc <- 2 * beta_mle - beta_b_mean
  
  # Predict with bias-corrected coefficients
  X <- model.matrix(formula, data = data)
  eta_bc <- X %*% beta_bc
  ps_bc <- 1 / (1 + exp(-eta_bc))
  ps_bc <- pmin(pmax(ps_bc, 1e-6), 1 - 1e-6)
  
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps_bc, (1 - p_marginal) / (1 - ps_bc))
  
  return(list(
    ps = as.vector(ps_bc),
    weights = weights,
    converged = TRUE,
    method = "BC-MLE"
  ))
}


#' Estimate propensity scores: Firth's penalized likelihood (FL)
#' 
#' Penalizes likelihood with Jeffreys prior; finite under separation.
#' Biases average predicted probability toward 0.5 (problematic for weighting).
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_firth <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  fit <- tryCatch(
    brglm2::glm(formula, data = data, family = binomial(), method = "brglmFit"),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "FL"))
  }
  
  ps <- predict(fit, type = "response")
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = ps,
    weights = weights,
    converged = TRUE,
    method = "FL"
  ))
}


#' Estimate propensity scores: Firth with Intercept Correction (FLIC)
#' 
#' Firth slopes + re-estimated intercept to restore unbiased average predicted probability.
#' Reference: Puhr et al. (2017)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_flic <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  # Step 1: Fit Firth
  fit_firth <- tryCatch(
    brglm2::glm(formula, data = data, family = binomial(), method = "brglmFit"),
    error = function(e) NULL
  )
  
  if (is.null(fit_firth)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "FLIC"))
  }
  
  beta_firth <- coef(fit_firth)
  
  # Step 2: Form linear predictor WITHOUT intercept
  X <- model.matrix(formula, data = data)
  eta_no_int <- X[, -1, drop = FALSE] %*% beta_firth[-1]
  
  # Step 3: Re-estimate intercept by MLE with eta_no_int as offset
  data$offset_flic <- eta_no_int
  fit_int <- tryCatch(
    glm(A_base ~ 1 + offset(offset_flic), data = data, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit_int) || !fit_int$converged) {
    # Fall back to Firth if intercept re-estimation fails
    ps <- predict(fit_firth, type = "response")
  } else {
    gamma_0 <- coef(fit_int)[1]
    eta_flic <- gamma_0 + eta_no_int
    ps <- 1 / (1 + exp(-eta_flic))
  }
  
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = as.vector(ps),
    weights = weights,
    converged = TRUE,
    method = "FLIC"
  ))
}


#' Estimate propensity scores: Firth with Added Covariate (FLAC)
#' 
#' Data augmentation approach; unbiased average PS, scale-invariant, low RMSE.
#' Reference: Puhr et al. (2017)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_flac <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  # Step 1: Fit Firth to get hat matrix
  fit_firth <- tryCatch(
    brglm2::glm(formula, data = data, family = binomial(), method = "brglmFit"),
    error = function(e) NULL
  )
  
  if (is.null(fit_firth)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "FLAC"))
  }
  
  # Hat diagonals
  X <- model.matrix(formula, data = data)
  W <- diag(fitted(fit_firth) * (1 - fitted(fit_firth)))
  H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  h_i <- diag(H)
  
  # Step 2: Build augmented dataset
  # Block (i): original data, weight 1, g=0
  # Block (ii): original data, weight h_i/2, g=1
  # Block (iii): original data with outcome flipped, weight h_i/2, g=1
  
  data_aug <- bind_rows(
    data %>% mutate(weight_aug = 1, g = 0),
    data %>% mutate(weight_aug = h_i / 2, g = 1),
    data %>% mutate(A_base_flip = 1 - A_base, weight_aug = h_i / 2, g = 1) %>%
      select(-A_base) %>%
      rename(A_base = A_base_flip)
  )
  
  # Step 3: Fit MLE logistic regression including g as covariate
  formula_aug <- update(formula, ~ . + g)
  fit_flac <- tryCatch(
    glm(formula_aug, data = data_aug, family = binomial(), weights = weight_aug),
    error = function(e) NULL
  )
  
  if (is.null(fit_flac) || !fit_flac$converged) {
    ps <- predict(fit_firth, type = "response")
  } else {
    # Predict for ORIGINAL data (g=0)
    data_pred <- data %>% mutate(g = 0)
    ps <- predict(fit_flac, newdata = data_pred, type = "response")
  }
  
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = as.vector(ps),
    weights = weights,
    converged = TRUE,
    method = "FLAC"
  ))
}


#' Estimate propensity scores: log-F(1,1) penalization (LF)
#' 
#' Independent log-F(1,1) priors on non-intercept coefficients.
#' Interpretable default; omitting intercept from penalty yields unbiased average PS.
#' Reference: Greenland & Mansournia (2015)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_logf <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  fit <- tryCatch(
    logistf::logistf(formula, data = data, firth = TRUE, pl = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "LF"))
  }
  
  # Predict
  X <- model.matrix(formula, data = data)
  eta <- X %*% fit$coefficients
  ps <- 1 / (1 + exp(-eta))
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = as.vector(ps),
    weights = weights,
    converged = TRUE,
    method = "LF"
  ))
}


#' Estimate propensity scores: Cauchy prior (CP)
#' 
#' Weakly informative Cauchy priors via bayesglm.
#' Near-zero PS bias, valid CIs.
#' Reference: Gelman et al. (2008)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_cauchy <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  fit <- tryCatch(
    arm::bayesglm(formula, data = data, family = binomial(),
                  prior.scale = 2.5, prior.scale.for.intercept = 10),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "CP"))
  }
  
  ps <- predict(fit, type = "response")
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = ps,
    weights = weights,
    converged = TRUE,
    method = "CP"
  ))
}


#' Estimate propensity scores: LASSO (L1 penalization)
#' 
#' L1 penalty with cross-validated lambda; variable selection; finite under separation.
#' Reference: Tibshirani (1996)
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return List with ps, weights, converged, method
#' @export
ps_lasso <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  # Prepare design matrix (excluding intercept, glmnet adds it)
  X <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  y <- data$A_base
  
  fit <- tryCatch(
    glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(ps = rep(NA, nrow(data)), weights = rep(NA, nrow(data)),
                converged = FALSE, method = "LASSO"))
  }
  
  # Predict at lambda.min
  ps <- as.vector(predict(fit, newx = X, s = "lambda.min", type = "response"))
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  p_marginal <- mean(data$A_base, na.rm = TRUE)
  weights <- ifelse(data$A_base == 1, p_marginal / ps, (1 - p_marginal) / (1 - ps))
  
  return(list(
    ps = ps,
    weights = weights,
    converged = TRUE,
    method = "LASSO"
  ))
}


#' Estimate all 8 PS methods for a single dataset
#' 
#' Convenience wrapper to run all estimators and return a named list.
#' 
#' @param data Person-trial dataset
#' @param formula Formula for PS model
#' @return Named list of results from each estimator
#' @export
estimate_all_ps <- function(data, formula = A_base ~ Lc_base + Lb_base + V) {
  
  list(
    MLE     = ps_mle(data, formula),
    BC_MLE  = ps_bc_mle(data, formula, B = 200),
    FL      = ps_firth(data, formula),
    FLIC    = ps_flic(data, formula),
    FLAC    = ps_flac(data, formula),
    LF      = ps_logf(data, formula),
    CP      = ps_cauchy(data, formula),
    LASSO   = ps_lasso(data, formula)
  )
}
