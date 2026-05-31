ps_model_wt <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c("sexM", "age_lt30", "age_30_39",
                                         "age_40_49", "age_ge50", "rural",
                                         "year", "po_preOUD", "phy_att_gr2",
                                         "social_assist", "ret_length", "home_grp_b",
                                         "incar_grp_bl", "drug_acute", "cci_gt0",
                                         "copd_asthma", "sed_psy_med", "po", "sud",
                                         "aud", "smd", "hcv", "tob", "udt", "ddose_prev",
                                         "clin_size", "virtual_care", "pss"),
                        continuous_vars = c("ddose_prev", "year", "ret_length")) {

  library(dplyr)
  library(nnet)

  model_data <- dat %>%
    dplyr::select(all_of(c(out_vars, clust_vars, exposure, baselinevars))) %>%
    mutate(across(all_of(setdiff(baselinevars, continuous_vars)), as.factor)) %>%
    filter(if_all(all_of(c(exposure, baselinevars)), ~ !is.na(.)))

  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))

  # Ensure exposure is a factor
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])

  # Fit multinomial PS model
  ps_fit <- multinom(formula_ps, data = model_data, trace = FALSE)
  ps_probs <- predict(ps_fit, type = "probs")

  ps_probs <- matrix(ps_probs, ncol = length(levels_exp))
  colnames(ps_probs) <- levels_exp

  # Get marginal probabilities (stabilizer)
  marginal_probs <- prop.table(table(model_data[[exposure]]))

  # Compute weights
  obs_exposure <- model_data[[exposure]]
  model_data$ps <- mapply(function(i, trt) ps_probs[i, trt],
                          i = seq_along(obs_exposure),
                          trt = obs_exposure)
  model_data$stab_numerator <- marginal_probs[as.character(obs_exposure)]
  model_data$wt_stab <- as.numeric(as.character(model_data$stab_numerator / model_data$ps))

  # Truncate weights
  p_low <- quantile(model_data$wt_stab, trunc_low)
  p_high <- quantile(model_data$wt_stab, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt_stab, p_low), p_high)

  # Truncate weights
  ps_probs_trunc <- as.data.frame(
    lapply(ps_probs, function(x) {
      p_low  <- quantile(x, trunc_low)
      p_high <- quantile(x, trunc_high)
      pmin(pmax(x, p_low), p_high)
    })
  )

  p_low <- quantile(model_data$ps, trunc_low)
  p_high <- quantile(model_data$ps, trunc_high)
  model_data$ps_trunc <- pmin(pmax(model_data$ps, p_low), p_high)

  ## ow
  tilting_fn <- apply(ps_probs, 1, function(x) sum(1/x))

  tilting_fn_trunc <- apply(ps_probs_trunc, 1, function(x) sum(1/x))

  model_data$ow_wt <- as.numeric(1/ (tilting_fn*model_data$ps))

  model_data$ow_wt_trunc <- as.numeric(1/ (tilting_fn_trunc*model_data$ps_trunc))

  return(list(model = ps_fit, data_weights = model_data))
}

ps_model_lasso <- function(dat, exposure = "trt",
                           trunc_low  = 0.01, trunc_high = 0.99,
                           out_vars, clust_vars,
                           baselinevars,
                           continuous_vars = c("ddose_prev", "year", "ret_length"),
                           alpha = 1,      # 1 = LASSO, 0 = ridge, 0.5 = elastic net
                           nfolds = 5) {   # cross-validation folds
  
  library(dplyr)
  library(glmnet)
  
  # Prepare data
  model_data <- dat %>%
    dplyr::select(all_of(c(out_vars, clust_vars, exposure, baselinevars))) %>%
    mutate(across(all_of(setdiff(baselinevars, continuous_vars)), as.factor)) %>%
    filter(if_all(all_of(c(exposure, baselinevars)), ~ !is.na(.)))
  
  valid_baselinevars <- baselinevars[
    sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)
  ]
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  
  # Prepare design matrix (glmnet requires matrix input)
  X <- model.matrix(formula_ps, data = model_data)[, -1]  # remove intercept
  y <- model_data[[exposure]]
  
  cat(sprintf("Fitting LASSO multinomial (alpha=%.1f, %d-fold CV)...\n", alpha, nfolds))
  
  # Fit with cross-validation
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "multinomial",
    alpha = alpha,
    nfolds = nfolds,
    type.measure = "deviance"
  )
  
  # Predict at lambda.min (or lambda.1se for more shrinkage)
  ps_probs <- predict(cv_fit, newx = X, s = "lambda.min", type = "response")
  ps_probs <- ps_probs[, , 1]  # Extract the probability matrix
  
  # Ensure column names match exposure levels
  if (ncol(ps_probs) != length(levels_exp)) {
    stop("LASSO output dimensions don't match exposure levels.")
  }
  colnames(ps_probs) <- levels_exp
  
  # Bound probabilities
  ps_probs <- apply(ps_probs, 2, function(x) pmin(pmax(x, 1e-6), 1 - 1e-6))
  ps_probs <- ps_probs / rowSums(ps_probs)  # renormalize
  
  # Compute weights (same as ps_model_wt)
  marginal_probs <- prop.table(table(model_data[[exposure]]))
  obs_exposure   <- model_data[[exposure]]
  
  model_data$ps <- mapply(function(i, trt) ps_probs[i, trt],
                          i = seq_along(obs_exposure), trt = obs_exposure)
  model_data$stab_numerator <- marginal_probs[as.character(obs_exposure)]
  model_data$wt_stab        <- as.numeric(model_data$stab_numerator / model_data$ps)
  
  # Truncated stabilised IPTW
  p_low  <- quantile(model_data$wt_stab, trunc_low)
  p_high <- quantile(model_data$wt_stab, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt_stab, p_low), p_high)
  
  # Truncated PS
  p_low  <- quantile(model_data$ps, trunc_low)
  p_high <- quantile(model_data$ps, trunc_high)
  model_data$ps_trunc <- pmin(pmax(model_data$ps, p_low), p_high)
  
  # Truncated PS probabilities (for OW)
  ps_probs_trunc <- apply(ps_probs, 2, function(x)
    pmin(pmax(x, quantile(x, trunc_low)), quantile(x, trunc_high)))
  
  # Overlap weights
  tilting_fn       <- apply(ps_probs,       1, function(x) sum(1 / x))
  tilting_fn_trunc <- apply(ps_probs_trunc, 1, function(x) sum(1 / x))
  
  model_data$ow_wt       <- as.numeric(1 / (tilting_fn       * model_data$ps))
  model_data$ow_wt_trunc <- as.numeric(1 / (tilting_fn_trunc * model_data$ps_trunc))
  
  return(list(
    model = cv_fit,
    data_weights = model_data,
    lambda_min = cv_fit$lambda.min
  ))
}


ps_rds = paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds")
names(dat)
dat %>% count(trt)

make_table1 <- function(ps_rds, bvars, cont_vars, clustered = FALSE) {
  
  dat <- readRDS(ps_rds)$data_weights |>
    mutate(incar_grp_bl = if_else(incar_grp_bl == "1_pastyear", 1L, 0L))
  
  factor_vars <- setdiff(bvars, cont_vars)
  
  # Unweighted
  tbl_before <- CreateTableOne(
    vars       = bvars,
    factorVars = factor_vars,
    strata     = "trt",
    data       = dat,
    addOverall = TRUE
  ) |>
    print(smd = TRUE, showAllLevels = FALSE, nonnormal = cont_vars,
          quote = FALSE, noSpaces = TRUE, printToggle = FALSE, test = FALSE)
  
  # Weighted (survey design)
  design <- if (clustered) {
    svydesign(id = ~studyid_rev, weights = ~dat$wt_trunc, data = dat)
  } else {
    svydesign(id = ~1, weights = ~dat$wt_trunc, data = dat)
  }
  
  tbl_after <- svyCreateTableOne(
    vars       = bvars,
    factorVars = factor_vars,
    strata     = "trt",
    data       = design,
    addOverall = FALSE
  ) |>
    print(smd = TRUE, showAllLevels = FALSE, nonnormal = cont_vars,
          quote = FALSE, noSpaces = TRUE, printToggle = FALSE, test = FALSE)
  
  cbind(tbl_before, SMD_weighted = tbl_after[, "SMD"])
}


# Helper: build one love plot panel for a single pairwise comparison
make_love_panel <- function(dat, bvars, var_labels, title_txt,
                            show_y_labels = TRUE) {
  
  fmla <- as.formula(paste("trt ~", paste(bvars, collapse = " + ")))
  
  prep <- bal.tab(
    fmla,
    data        = dat,
    weights     = dat$wt_trunc,
    estimand    = "ATE",
    m.threshold = 0.1,
    s.d.denom   = "pooled",
    un          = TRUE,
    stats       = "mean.diffs"
  )
  
  # Alternating row background (matches number of covariates in bal.tab output)
  n_vars  <- nrow(prep$Balance)
  bg_data <- data.frame(
    var_index  = seq_len(n_vars),
    fill_color = rep(c("gray65", "white"), length.out = n_vars)
  )
  
  plt <- love.plot(
    prep,
    stat         = "mean.diffs",
    stars        = "none",
    thresholds   = c(m = 0.1),
    abs          = FALSE,
    var.names    = var_labels,
    sample.names = c("Before weighting", "After weighting"),
    position     = "none",
    line         = FALSE,
    colors       = c("black", "red"),
    shapes       = c(16, 16),
    title        = title_txt
  ) +
    geom_rect(
      data        = bg_data,
      aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
          xmin = -Inf, xmax = Inf, fill = fill_color),
      inherit.aes = FALSE,
      alpha       = 0.2
    ) +
    scale_x_continuous(limits = c(x_low, x_high),
                       breaks = x_breaks, expand = c(0, 0)) +
    scale_fill_identity() +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(x = "SMD", color = "", shape = "") +
    theme(
      title              = element_text(color = "black", size = 14),
      panel.grid         = element_blank(),
      axis.line          = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text.x        = element_text(color = "black", size = 13),
      axis.ticks.y       = element_line(linewidth = 1.2, color = "black"),
      axis.text.y        = if (show_y_labels)
        element_text(color = "black", size = 13)
      else
        element_blank(),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 13, color = "black")
    )
  
  plt
}

# Helper: build one 2-panel row (reduction | increase) for one trial × cohort
make_love_row <- function(ps_rds, bvars, var_labels, clustered = FALSE,
                          trt_recode = c("no reversion" = "No dose change",
                                         "reduction"    = "Dose reduction",
                                         "increase"     = "Dose increase")) {
  
  dat <- readRDS(ps_rds)$data_weights |>
    mutate(
      trt          = recode(as.character(trt), !!!trt_recode),
      incar_grp_bl = if_else(incar_grp_bl == "1_pastyear", 1L, 0L)
    )
  
  # Reduction vs no change
  d_red <- subset(dat, trt %in% c("No dose change", "Dose reduction"))
  d_red$trt <- droplevels(factor(d_red$trt))
  
  p_red <- make_love_panel(d_red, bvars, var_labels,
                           title_txt     = "Dose reduction\nversus no change",
                           show_y_labels = TRUE)
  
  # Increase vs no change
  d_inc <- subset(dat, trt %in% c("No dose change", "Dose increase"))
  d_inc$trt <- droplevels(factor(d_inc$trt))
  
  p_inc <- make_love_panel(d_inc, bvars, var_labels,
                           title_txt     = "Dose increase\nversus no change",
                           show_y_labels = FALSE)
  
  p_red + plot_spacer() + p_inc +
    plot_layout(widths = c(1, 0.07, 1), guides = "collect") &
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 13))
}



make_love_row_v1 <- function(ps_rds, bvars, var_labels, clustered = FALSE,
                          trt_recode = c("no reversion" = "No dose change",
                                         "1-<50%"    = "1-<50% reduction",
                                         ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                                         "increase"     = "Dose increase")) {
  
  dat <- readRDS(ps_rds)$data_weights |>
    mutate(
      trt          = recode(as.character(trt), !!!trt_recode),
      incar_grp_bl = if_else(incar_grp_bl == "1_pastyear", 1L, 0L)
    )
  
  # 1-<50% vs no change
  d_red_lt50 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
  d_red_lt50$trt <- droplevels(factor(d_red_lt50$trt))
  
  p_red_lt50 <- make_love_panel(d_red_lt50, bvars, var_labels,
                           title_txt     = "1-<50% reduction\nversus no change",
                           show_y_labels = TRUE)
  
  # >=50% reduction/initial dose* vs no change
  d_red_ge50 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
  d_red_ge50$trt <- droplevels(factor(d_red_ge50$trt))
  
  p_red_ge50 <- make_love_panel(d_red_ge50, bvars, var_labels,
                           title_txt     = ">=50% reduction/initial dose*\nversus no change",
                           show_y_labels = TRUE)
  
  
  
  # Increase vs no change
  d_inc <- subset(dat, trt %in% c("No dose change", "Dose increase"))
  d_inc$trt <- droplevels(factor(d_inc$trt))
  
  p_inc <- make_love_panel(d_inc, bvars, var_labels,
                           title_txt     = "Dose increase\nversus no change",
                           show_y_labels = FALSE)
  
  p_red_lt50 + plot_spacer()  + p_red_ge50 + plot_spacer() + p_inc +
    plot_layout(widths = c(1, 0.07, 1, 0.07, 1), guides = "collect") &
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 13))
}


expand_person_period <- function(Dat, time_var, event_var, interval_length = 7){
  Dat <- Dat %>%
    mutate(max_int = ceiling(.data[[time_var]] / interval_length))
  
  Dat_expanded <- Dat %>%
    tidyr::uncount(weights = max_int, .id = "interval", .remove = FALSE) %>%
    mutate(
      start = (interval - 1) * interval_length,
      stop  = interval * interval_length,
      event = ifelse(interval == max_int, .data[[event_var]], 0)
    )
  return(Dat_expanded)
}



predict_cuminc <- function(Dat_expanded, trt_var, covs, horizons,
                              cluster_var, samples = 1000, conf_int = TRUE,
                              ref_group = "no reversion", wt_type="wt_trunc") {
  
  # Formula for pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  # Fit GLM
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = Dat_expanded[[wt_type]])
  
  # Cluster-robust covariance
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  
  coefs <- coef(model)
  
  # Function to calculate risk at a given horizon
  risk_fun <- function(coefs, horizon) {
    nd <- Dat_expanded %>% filter(stop <= horizon)
    
    risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level) {
      nd_trt <- nd
      nd_trt[[trt_var]] <- trt_level
      
      mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
      p_trt <- plogis(mm_trt %*% coefs)
      nd_trt$p_event <- p_trt
      
      df_pred <- nd_trt %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = trt_level, horizon = horizon, risk = mean(df_pred$risk))
    })
    
    bind_rows(risks)
  }
  
  # Point estimates
  point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  if(!conf_int) return(point)
  
  # Simulate coefficient draws for CI
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    bind_rows(lapply(horizons, function(h) risk_fun(sim_coefs[i,], h)))
  })
  sim_df <- bind_rows(sim_list, .id = "sim")
  
  # Risk CI
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk),  # point estimate from simulation
      .groups = "drop"
    )
  
  # Reference group
  if(is.null(ref_group)) {
    ref_group <- levels(factor(Dat_expanded[[trt_var]]))[1]
  }
  
  # Compute RD and RR in each simulation
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(
      RD = risk - ref_risk,
      RR = risk / ref_risk
    )
  
  # Summarise RD and RR
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025),
      RD_upper = quantile(RD, 0.975),
      RD = mean(RD),
      RR_lower = quantile(RR, 0.025),
      RR_upper = quantile(RR, 0.975),
      RR = mean(RR),
      .groups = "drop"
    )
  
  list(
    risks = risk_ci,
    contrasts = contrast_ci,
    ref_group = ref_group
  )
}

