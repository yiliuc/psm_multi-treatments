###############################################################################
# (A)  Your existing helper functions for data generation, PS estimation, 
#      common support, KMC, and the various matching procedures:

simulate_multitreatment_data <- function(n = 3000, seed = 123) {
  # set.seed(seed)
  
  # Generate covariates
  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  X3 <- rnorm(n, mean = 0, sd = 1)
  X4 <- rnorm(n, mean = 0, sd = 1)
  
  # Define linear predictors for the multinomial logistic
  f1 <-  0.5 * X1 + 0.3 * X2 + 0.4 * X3 - 0.4 * X1^2
  f2 <- -0.2 * X1 + 0.6 * X2 - 0.3 * X3 + 0.3 * X2^2
  
  # Multinomial probabilities
  denom <- 1 + exp(f1) + exp(f2)
  p1 <- 1       / denom
  p2 <- exp(f1) / denom
  p3 <- exp(f2) / denom
  
  # Draw treatment T
  T <- numeric(n)
  for (i in seq_len(n)) {
    T[i] <- sample(x = 1:3, size = 1, prob = c(p1[i], p2[i], p3[i]))
  }
  
  # Outcome model + noise
  base_mean <- 4 + 0.3 * X1 - 0.2 * X2 + 0.6 * X4
  Y <- base_mean +
    ifelse(T == 2, -0.5, 0) +
    ifelse(T == 3,  0.7, 0) +
    rnorm(n, mean = 0, sd = 1)
  
  data.frame(
    id        = 1:n,
    X1        = X1,
    X2        = X2,
    X3        = X3,
    X4        = X4,
    Treatment = T,
    Y         = Y
  )
}
################################################################################################
library(nnet)  # for multinom
estimate_propensity_scores <- function(data) {
  data$Treatment <- as.factor(data$Treatment)
  model <- multinom(Treatment ~ . - Treatment - id - Y, data = data, trace = FALSE)
  propensity_scores <- predict(model, type = "probs")
  treatment_levels <- levels(data$Treatment)
  colnames(propensity_scores) <- paste0("PS", treatment_levels)
  return(propensity_scores)
}
################################################################################################
common_support <- function(ps_matrix, treatment_col) {
  treatment_groups <- unique(ps_matrix[[treatment_col]])
  low_values  <- numeric(length(treatment_groups))
  high_values <- numeric(length(treatment_groups))
  
  for (target_treatment in treatment_groups) {
    ps_target <- ps_matrix[[paste0("PS", target_treatment)]]
    
    group_min <- sapply(treatment_groups, function(actual_treatment) {
      group_data <- ps_matrix[ps_matrix[[treatment_col]] == actual_treatment, ]
      ps_values <- group_data[[paste0("PS", target_treatment)]]
      return(min(ps_values))
    })
    group_max <- sapply(treatment_groups, function(actual_treatment) {
      group_data <- ps_matrix[ps_matrix[[treatment_col]] == actual_treatment, ]
      ps_values <- group_data[[paste0("PS", target_treatment)]]
      return(max(ps_values))
    })
    
    # r(t,X)^low, r(t,X)^high
    low_values[target_treatment]  <- max(group_min)
    high_values[target_treatment] <- min(group_max)
  }
  
  list(
    Low  = low_values,
    High = high_values
  )
}

filter_valid_subjects_indices <- function(ps_data, css) {
  low  <- css$Low
  high <- css$High
  
  is_valid <- apply(ps_data, 1, function(row) {
    all(row >= low & row <= high)
  })
  which(is_valid)
}
################################################################################################
# Implementing KMC
perform_kmc <- function(ps_data, ref_treatment, t_prime, K = 5) {
  # Perform the KMC based on the vector of GPS exclude PS(t) and PS(t_prime)
  
  # Exclude the propensity scores for ref_treatment and t_prime
  ps_subset <- ps_data[, !colnames(ps_data) %in% c(ref_treatment, t_prime), drop = FALSE]
  
  # Check if there are enough columns to proceed
  if (ncol(ps_subset) < 1) {
    stop("Insufficient columns after excluding reference and target treatments. Ensure the data contains enough treatments.")
  }
  
  # Apply logit transformation to the remaining propensity scores
  logit_transform <- function(p) log(p / (1 - p))
  ps_logit <- apply(ps_subset, 2, logit_transform)
  
  # Ensure K does not exceed the number of rows in the data
  if (K > nrow(ps_logit)) {
    stop("Number of clusters (K) cannot exceed the number of rows in the data.")
  }
  
  # Perform KMC with the specified K
  set.seed(123)
  final_kmc <- kmeans(ps_logit, centers = K, nstart = 20)
  
  # Return cluster assignments and KMC model
  return(list(
    clusters = final_kmc$cluster,
    kmc_model = final_kmc
  ))
}
################################################################################################
library(MatchIt)
library(Matching)

perform_matching <- function(data, ref_treatment, t_prime, epsilon = 0.2) {
  all_matches <- data.frame()
  unique_clusters <- unique(data$Cluster)
  
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    # logit(PS(ref_treatment))
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    
    # Keep T in {ref_treatment, t_prime}
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    caliper_width <- epsilon * sd(cluster_data$logit_r_t, na.rm = TRUE)
    X  <- cluster_data$logit_r_t
    Tr <- cluster_data$treatment_indicator
    
    match_result <- Match(
      Y = NULL,
      Tr = Tr,
      X = X,
      M = 1,
      caliper = caliper_width,
      replace = TRUE
    )
    
    matched_id <- data.frame(
      id_ref     = cluster_data$id[ match_result$index.treated ],
      t_ref      = cluster_data$Treatment[ match_result$index.treated ],
      id_t_prime = cluster_data$id[ match_result$index.control ],
      t_prime    = cluster_data$Treatment[ match_result$index.control ]
    )
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################################
perform_matching_2 <- function(data, ref_treatment, t_prime, epsilon = 0.2, n_match = 2) {
  all_matches <- data.frame()
  unique_clusters <- unique(data$Cluster)
  
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    caliper_width <- epsilon * sd(cluster_data$logit_r_t, na.rm = TRUE)
    X  <- cluster_data$logit_r_t
    Tr <- cluster_data$treatment_indicator
    
    match_result <- Match(
      Y = NULL,
      Tr = Tr,
      X = X,
      M = n_match,
      caliper = caliper_width,
      replace = TRUE
    )
    
    matched_id <- data.frame(
      id_ref     = cluster_data$id[ match_result$index.treated ],
      t_ref      = cluster_data$Treatment[ match_result$index.treated ],
      id_t_prime = cluster_data$id[ match_result$index.control ],
      t_prime    = cluster_data$Treatment[ match_result$index.control ]
    )
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################################
perform_matching_vm_md <- function(data, ref_treatment, t_prime, epsilon = 0.2) {
  all_matches <- data.frame()
  unique_clusters <- unique(data$Cluster)
  
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t      <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    cluster_data$logit_r_tprime <- logit_transform(cluster_data[[paste0("PS", t_prime)]])
    
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    mat_2d <- as.matrix(cluster_data[, c("logit_r_t", "logit_r_tprime")])
    
    Sigma <- cov(mat_2d)
    Sigma_inv <- solve(Sigma)
    L <- chol(Sigma_inv)
    WeightedX <- t(L %*% t(mat_2d))
    
    # caliper in WeightedX space
    WeightedX_mean <- apply(WeightedX, 1, function(z) sqrt(sum(z^2)))
    caliper_width <- epsilon * sd(WeightedX_mean, na.rm = TRUE)
    
    match_result <- Match(
      Y = NULL,
      Tr = cluster_data$treatment_indicator,
      X = WeightedX,  
      M = 1,
      caliper = caliper_width,
      replace = TRUE
    )
    
    matched_id <- data.frame(
      id_ref     = cluster_data$id[ match_result$index.treated ],
      t_ref      = cluster_data$Treatment[ match_result$index.treated ],
      id_t_prime = cluster_data$id[ match_result$index.control ],
      t_prime    = cluster_data$Treatment[ match_result$index.control ],
      cluster    = k
    )
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################################
perform_matching_vm_mdnc <- function(data, ref_treatment, t_prime) {
  all_matches <- data.frame()
  unique_clusters <- unique(data$Cluster)
  
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t      <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    cluster_data$logit_r_tprime <- logit_transform(cluster_data[[paste0("PS", t_prime)]])
    
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    mat_2d <- as.matrix(cluster_data[, c("logit_r_t", "logit_r_tprime")])
    
    Sigma <- cov(mat_2d)
    Sigma_inv <- solve(Sigma)
    L <- chol(Sigma_inv)
    WeightedX <- t(L %*% t(mat_2d))
    
    # "No caliper" => caliper=0
    match_result <- Match(
      Y = NULL,
      Tr = cluster_data$treatment_indicator,
      X = WeightedX,
      M = 1,
      caliper = 0,
      replace = TRUE
    )
    
    matched_id <- data.frame(
      id_ref     = cluster_data$id[ match_result$index.treated ],
      t_ref      = cluster_data$Treatment[ match_result$index.treated ],
      id_t_prime = cluster_data$id[ match_result$index.control ],
      t_prime    = cluster_data$Treatment[ match_result$index.control ],
      cluster    = k
    )
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}

###############################################################################
# (B)  Function to compute standardized bias:

calc_standardized_bias <- function(data, covariates, ref_treatment = 1) {
  # This function expects data with columns:
  #   - id
  #   - Treatment
  #   - pair (or matched set index)
  #   - the covariates in 'covariates'
  #
  # Weighted approach: each subject i has frequency psi_i (# times matched).
  
  n_pairs <- length(unique(data$pair))
  count_by_id <- table(data$id)
  psi_lookup  <- setNames(as.numeric(count_by_id), names(count_by_id))
  get_psi <- function(id_i) psi_lookup[as.character(id_i)]
  
  trt_levels <- sort(unique(data$Treatment))
  if (length(trt_levels) != 2) {
    stop("Data must have exactly 2 treatments for this SB calculation!")
  }
  
  # Weighted means
  wmean_mat <- matrix(
    NA, 
    nrow = length(covariates), 
    ncol = 2, 
    dimnames = list(covariates, paste0("T", trt_levels))
  )
  
  for (p in seq_along(covariates)) {
    cov_name <- covariates[p]
    for (trt in trt_levels) {
      rows_trt <- which(data$Treatment == trt)
      if (length(rows_trt) == 0) {
        wmean_mat[p, paste0("T", trt)] <- NA
        next
      }
      Xp_vals  <- data[[cov_name]][rows_trt]
      psi_vals <- sapply(data$id[rows_trt], get_psi)
      weighted_sum <- sum(Xp_vals * psi_vals)
      wmean_mat[p, paste0("T", trt)] <- weighted_sum / n_pairs
    }
  }
  
  # Weighted SD in the reference group
  wsd_vec <- numeric(length(covariates))
  names(wsd_vec) <- covariates
  
  ref_col <- paste0("T", ref_treatment)
  rows_ref <- which(data$Treatment == ref_treatment)
  psi_vals_ref <- sapply(data$id[rows_ref], get_psi)
  
  for (p in seq_along(covariates)) {
    Xp_vals_ref <- data[[covariates[p]]][rows_ref]
    wmean_ref   <- wmean_mat[p, ref_col]
    
    var_num <- sum(psi_vals_ref * (Xp_vals_ref - wmean_ref)^2)
    var_den <- sum(psi_vals_ref)
    if (var_den > 1e-10) {
      wsd_vec[p] <- sqrt(var_num / var_den)
    } else {
      wsd_vec[p] <- NA
    }
  }
  
  # Standardized difference for the other group vs. reference
  other_treatment <- setdiff(trt_levels, ref_treatment)
  other_col <- paste0("T", other_treatment)
  
  sb_vec <- numeric(length(covariates))
  names(sb_vec) <- covariates
  for (p in seq_along(covariates)) {
    denom <- wsd_vec[p]
    if (is.na(denom) || denom < 1e-10) {
      sb_vec[p] <- NA
    } else {
      mean_other <- wmean_mat[p, other_col]
      mean_ref   <- wmean_mat[p, ref_col]
      sb_vec[p]  <- (mean_other - mean_ref) / denom
    }
  }
  
  abs_sb <- abs(sb_vec)
  return(abs_sb)
}
###############################################################################
calc_ate_ipw <- function(data, t1 = 1, t2 = 2,
                         y_col = "Y",
                         trt_col = "Treatment",
                         ps_cols = c("PS1", "PS2", "PS3")) {
  # data must contain:
  #  - Y (outcome)
  #  - Treatment in {1,2,3} (the actual treatment each row got)
  #  - PS1, PS2, PS3 = propensity scores for T=1, T=2, T=3.
  #
  # t1, t2: which two treatments do you want to compare?
  
  # 1) Compute inverse-propensity weights for each row
  #    If row i got T=1, weight = 1/PS1[i].
  #    If row i got T=2, weight = 1/PS2[i].
  #    If row i got T=3, weight = 1/PS3[i].
  
  data$ipw <- apply(data, 1, function(row) {
    trt_val <- as.numeric(row[[trt_col]])
    if (trt_val == 1) {
      return(1 / as.numeric(row[[ ps_cols[1] ]]))  # 1/PS1
    } else if (trt_val == 2) {
      return(1 / as.numeric(row[[ ps_cols[2] ]]))  # 1/PS2
    } else if (trt_val == 3) {
      return(1 / as.numeric(row[[ ps_cols[3] ]]))  # 1/PS3
    } else {
      return(NA)
    }
  })
  
  # 2) For each of t1, t2, compute the weighted mean of Y 
  #    using the subset that actually received t1 or t2.
  
  # Weighted mean function
  wmean <- function(x, w) sum(x * w) / sum(w)
  
  # Subset for t1
  idx_t1 <- data[[trt_col]] == t1
  mu_t1 <- wmean(
    x = data[[y_col]][idx_t1],
    w = data$ipw[idx_t1]
  )
  
  # Subset for t2
  idx_t2 <- data[[trt_col]] == t2
  mu_t2 <- wmean(
    x = data[[y_col]][idx_t2],
    w = data$ipw[idx_t2]
  )
  
  # 3) Estimated ATE = mu_t1 - mu_t2
  ate_t1_t2 <- mu_t1 - mu_t2
  
  # Return a small summary in a list
  return(list(
    mu_t1      = mu_t1,
    mu_t2      = mu_t2,
    ATE_t1_t2  = ate_t1_t2
  ))
}
###############################################################################
# (C)  A single simulation that returns the final measure \bar{Max2SB} 
#      for one matching method:

library(dplyr)

one_simulation <- function(
    n                = 3000,
    K                = 3,
    matching_method  = c("M1", "M2", "VM_MD"),
    epsilon          = 0.2,   
    n_match          = 2,     # used only for "M2"
    seed             = NULL,
    covariates_for_SB = c("X1","X2","X3","X4")
) {
  if(!is.null(seed)) set.seed(seed)
  
  ## 1) Simulate data
  data <- simulate_multitreatment_data(n = n)
  
  ## 2) Estimate PS & drop out-of-support
  ps <- estimate_propensity_scores(data)
  ps_trt <- as.data.frame(cbind(ps, Treatment = data$Treatment))
  css <- common_support(ps_trt, "Treatment")
  valid_index <- filter_valid_subjects_indices(ps, css)
  valid_data  <- data[valid_index,]
  
  ps_new <- estimate_propensity_scores(valid_data)
  valid_data_ps <- cbind(valid_data, ps_new)
  
  ## 3) For each pair of treatments (1 vs 2), (1 vs 3), (2 vs 3):
  ##    (a) KMC => cluster assignments
  ##    (b) call the chosen matching function => matched pairs
  ##    (c) build final matched cohort => compute standardized biases
  
  # Helper: choose the correct matching function
  run_matching_func <- function(method, data_input, ref_t, t_prime_t) {
    if (method == "M1") {
      out <- perform_matching(
        data          = data_input, 
        ref_treatment = ref_t,
        t_prime       = t_prime_t,
        epsilon       = epsilon
      )
    } else if (method == "M2") {
      out <- perform_matching_2(
        data          = data_input,
        ref_treatment = ref_t,
        t_prime       = t_prime_t,
        epsilon       = epsilon,
        n_match       = n_match
      )
    } else if (method == "VM_MD") {
      out <- perform_matching_vm_md(
        data          = data_input,
        ref_treatment = ref_t,
        t_prime       = t_prime_t,
        epsilon       = epsilon
      )
    } else if (method == "VM_MDnc") {
      out <- perform_matching_vm_mdnc(
        data          = data_input,
        ref_treatment = ref_t,
        t_prime       = t_prime_t
      )
    } else {
      stop("Unrecognized matching method!")
    }
    return(out)
  }
  
  ###########
  # Pair (1 vs 2)
  kmc_12 <- perform_kmc(
    ps_data       = ps_new,
    ref_treatment = "PS1",
    t_prime       = "PS2",
    K             = K
  )
  valid_data_ps_12 <- cbind(valid_data_ps, Cluster = kmc_12$clusters)
  matches_1_2 <- run_matching_func(matching_method, valid_data_ps_12, ref_t=1, t_prime=2)
  
  final_cohort_1_2 <- data.frame(
    id   = c(rbind(matches_1_2$id_ref, matches_1_2$id_t_prime)),
    t    = c(rbind(matches_1_2$t_ref,  matches_1_2$t_prime )),
    pair = rep(seq_len(nrow(matches_1_2)), each = 2)
  )
  final_cohort_1_2 <- dplyr::left_join(
    final_cohort_1_2,
    valid_data_ps,
    by = c("id", "t" = "Treatment")
  ) %>%
    dplyr::rename(Treatment = t)
  
  abs_sb_12 <- calc_standardized_bias(
    data         = final_cohort_1_2,
    covariates   = covariates_for_SB,
    ref_treatment = 1
  )
  
  ###########
  # Pair (1 vs 3)
  kmc_13 <- perform_kmc(
    ps_data       = ps_new,
    ref_treatment = "PS1",
    t_prime       = "PS3",
    K             = K
  )
  valid_data_ps_13 <- cbind(valid_data_ps, Cluster = kmc_13$clusters)
  matches_1_3 <- run_matching_func(matching_method, valid_data_ps_13, ref_t=1, t_prime=3)
  
  final_cohort_1_3 <- data.frame(
    id   = c(rbind(matches_1_3$id_ref, matches_1_3$id_t_prime)),
    t    = c(rbind(matches_1_3$t_ref,  matches_1_3$t_prime )),
    pair = rep(seq_len(nrow(matches_1_3)), each = 2)
  )
  final_cohort_1_3 <- dplyr::left_join(
    final_cohort_1_3,
    valid_data_ps,
    by = c("id", "t" = "Treatment")
  ) %>%
    dplyr::rename(Treatment = t)
  
  abs_sb_13 <- calc_standardized_bias(
    data         = final_cohort_1_3,
    covariates   = covariates_for_SB,
    ref_treatment = 1
  )
  
  ###########
  # Pair (2 vs 3)
  kmc_23 <- perform_kmc(
    ps_data       = ps_new,
    ref_treatment = "PS2",
    t_prime       = "PS3",
    K             = K
  )
  valid_data_ps_23 <- cbind(valid_data_ps, Cluster = kmc_23$clusters)
  matches_2_3 <- run_matching_func(matching_method, valid_data_ps_23, ref_t=2, t_prime=3)
  
  final_cohort_2_3 <- data.frame(
    id   = c(rbind(matches_2_3$id_ref, matches_2_3$id_t_prime)),
    t    = c(rbind(matches_2_3$t_ref,  matches_2_3$t_prime )),
    pair = rep(seq_len(nrow(matches_2_3)), each = 2)
  )
  final_cohort_2_3 <- dplyr::left_join(
    final_cohort_2_3,
    valid_data_ps,
    by = c("id", "t" = "Treatment")
  ) %>%
    dplyr::rename(Treatment = t)
  
  abs_sb_23 <- calc_standardized_bias(
    data         = final_cohort_2_3,
    covariates   = covariates_for_SB,
    ref_treatment = 2
  )
  
  ## 4) Combine the standardized biases
  stopifnot(all(names(abs_sb_12) == names(abs_sb_13)))
  stopifnot(all(names(abs_sb_12) == names(abs_sb_23)))
  
  covariate_names <- names(abs_sb_12)
  max2sb <- numeric(length(covariate_names))
  for(i in seq_along(covariate_names)) {
    p_name <- covariate_names[i]
    # max over the absolute SB for (1,2), (1,3), (2,3)
    max2sb[i] <- max(
      abs_sb_12[p_name],
      abs_sb_13[p_name],
      abs_sb_23[p_name],
      na.rm = TRUE
    )
  }
  mean_max2sb <- mean(max2sb, na.rm = TRUE)
  
  ## 5) Also compute ATE(1,2) and ATE(1,3) from the *matched data*
  #    using your provided IPW approach.
  
  ate_12_list <- calc_ate_ipw(
    data    = final_cohort_1_2,
    t1      = 1,
    t2      = 2,
    y_col   = "Y",
    trt_col = "Treatment",
    ps_cols = c("PS1","PS2","PS3")
  )
  ate_13_list <- calc_ate_ipw(
    data    = final_cohort_1_3,
    t1      = 1,
    t2      = 3,
    y_col   = "Y",
    trt_col = "Treatment",
    ps_cols = c("PS1","PS2","PS3")
  )
  
  # We'll return just the numerical ATE from each
  ate_12 <- ate_12_list$ATE_t1_t2
  ate_13 <- ate_13_list$ATE_t1_t2
  
  # 6) Compute matching rate:
  #    A subject is considered matched if it appears in either final_cohort_1_2 or final_cohort_1_3.
  #    The total sample size was n=3000, so:
  matched_ids_12 <- unique(final_cohort_1_2$id)
  matched_ids_13 <- unique(final_cohort_1_3$id)
  matched_ids_23 <- unique(final_cohort_2_3$id)
  matched_union  <- unique(c(matched_ids_12, matched_ids_13, matched_ids_23))
  matching_rate  <- length(matched_union) / n  # proportion matched
  
  ## Return everything in a list
  return(list(
    max2sb = mean_max2sb,
    ate_12 = ate_12,
    ate_13 = ate_13,
    matching_rate = matching_rate
  ))
}


###############################################################################
# (D)  Finally, run B simulations and collect the results for each matching method

run_simulations <- function(B, n = 3000, K = 3) {
  # For M1
  results_M1_max2sb  <- numeric(B)
  results_M1_ate12   <- numeric(B)
  results_M1_ate13   <- numeric(B)
  results_M1_matchrate <- numeric(B)
  
  # For M2
  results_M2_max2sb  <- numeric(B)
  results_M2_ate12   <- numeric(B)
  results_M2_ate13   <- numeric(B)
  results_M2_matchrate <- numeric(B)
  
  # For VM_MD
  results_VM_MD_max2sb  <- numeric(B)
  results_VM_MD_ate12   <- numeric(B)
  results_VM_MD_ate13   <- numeric(B)
  results_VM_MD_matchrate <- numeric(B)
  
  for (b in seq_len(B)) {
    seed_b <- 123 + b  # optional for reproducibility
    
    # --- M1 ---
    sim_M1 <- one_simulation(
      n = n, K = K,
      matching_method  = "M1",
      seed            = seed_b
    )
    results_M1_max2sb[b]   <- sim_M1$max2sb
    results_M1_ate12[b]    <- sim_M1$ate_12
    results_M1_ate13[b]    <- sim_M1$ate_13
    results_M1_matchrate[b]<- sim_M1$matching_rate
    
    # --- M2 ---
    sim_M2 <- one_simulation(
      n = n, K = K,
      matching_method  = "M2",
      seed            = seed_b
    )
    results_M2_max2sb[b]   <- sim_M2$max2sb
    results_M2_ate12[b]    <- sim_M2$ate_12
    results_M2_ate13[b]    <- sim_M2$ate_13
    results_M2_matchrate[b]<- sim_M2$matching_rate
    
    # --- VM_MD ---
    sim_VM <- one_simulation(
      n = n, K = K,
      matching_method  = "VM_MD",
      seed            = seed_b
    )
    results_VM_MD_max2sb[b]   <- sim_VM$max2sb
    results_VM_MD_ate12[b]    <- sim_VM$ate_12
    results_VM_MD_ate13[b]    <- sim_VM$ate_13
    results_VM_MD_matchrate[b]<- sim_VM$matching_rate
    
    if (b %% 20 == 0) {
      cat("Finished iteration:", b, "of", B, "\n")
    }
  }
  
  return(list(
    VM = list(
      max2sb     = results_M1_max2sb,
      ate_12     = results_M1_ate12,
      ate_13     = results_M1_ate13,
      match_rate = results_M1_matchrate
    ),
    VM2 = list(
      max2sb     = results_M2_max2sb,
      ate_12     = results_M2_ate12,
      ate_13     = results_M2_ate13,
      match_rate = results_M2_matchrate
    ),
    VM_MD = list(
      max2sb     = results_VM_MD_max2sb,
      ate_12     = results_VM_MD_ate12,
      ate_13     = results_VM_MD_ate13,
      match_rate = results_VM_MD_matchrate
    )
  ))
}
###############################################################################
# (E)  Example usage
system.time({
  sim_results <- run_simulations(B=300, n=3000, K=3)
})
save(sim_results, file = "sim300_4var_s2_results.Rdata")
# We then get four vectors in sim_results$M1, sim_results$M2, etc.
# each containing the distribution of \bar{Max2SB} over B runs.

# You can then summarize each vector, e.g.:
# mean(sim_results$M1)
# sd(sim_results$M1)
# etc.
###############################################################################
