################################################################################
simulate_multitreatment_data <- function(n = 1000, seed = 123) {
  set.seed(seed)
  
  # 1. Generate covariates:
  #    - X1, X2, X6 ~ N(0,1)
  #    - X3, X5, X7 ~ Bernoulli(0.5)
  #    - X4 ~ N(0,1)
  #    (You can adjust distributions as needed.)
  
  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  X3 <- rbinom(n, size = 1, prob = 0.5)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  X5 <- rnorm(n, mean = 0, sd = 1)
  X6 <- rbinom(n, size = 1, prob = 0.5)
  X7 <- rnorm(n, mean = 0, sd = 1)
  X8 <- rbinom(n, size = 1, prob = 0.5)
  
  # 2. Define linear predictors for the multinomial logistic:
  #    f1(X) and f2(X) in your formula.
  
  f1 <- 0.5*X1 + 0.3*X2 - 0.3*X3 + 0.6*X4 - 0.25*X5 + 0.3*X6
  f2 <- -0.3*X1 + 0.5*X2 + 0.2*X3 - 0.3*X4 - 0.5*X5 + 0.2*X6
  
  # 3. Compute multinomial probabilities for T = 1, 2, 3:
  
  denom <- 1 + exp(f1) + exp(f2)
  p1 <- 1 / denom
  p2 <- exp(f1) / denom
  p3 <- exp(f2) / denom
  
  # 4. Draw one treatment T (1, 2, or 3) per subject using those probabilities:
  
  T <- numeric(n)
  for (i in seq_len(n)) {
    T[i] <- sample(x = 1:3, size = 1, prob = c(p1[i], p2[i], p3[i]))
  }
  
  # 5. Define the outcome model:
  #    E[Y | X] = 4 + 0.3X1 - 0.4X2 - 0.7X3 + 0.6X6 - 0.3X7
  #               - 0.5 * I{T=2} + 0.7 * I{T=3}
  #    Add random Normal(0,1) noise to get Y.
  
  base_mean <- 4 + 0.3*X1 - 0.4*X2 - 0.7*X3 + 0.2 * X4 + 0.6*X7 - 0.3*X8
  Y <- base_mean +
    ifelse(T == 2, -0.5, 0) +
    ifelse(T == 3,  0.7, 0)
  
  # Return a data frame with everything
  return(data.frame(id = 1:n,
                    X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8,
                    Treatment = T,
                    Y = Y
  ))
}

# Example usage
# data <- simulate_multitreatment_data(3000)
# 
# t1 <- final_cohort[final_cohort$Treatment == 1,]
# t2 <- final_cohort[final_cohort$Treatment == 2,]
# t3 <- final_cohort[final_cohort$Treatment == 3,]
# mean(t2$Y-t1$Y)
# mean(t3$Y-t1$Y)
# mean(t3$Y-t2$Y)

# List of scenarios
# f1 <- 0.5*X1 + 0.3*X2 - 0.3*X3 + 0.6*X4 - 0.25*X5 + 0.3*X6 - 0.4*X1^2
# f2 <- -0.3*X1 + 0.5*X2 + 0.2*X3 - 0.3*X4 - 0.5*X5 + 0.2*X6 + 0.3*X2^2
# 
# f1 <- 0.5*X1 + 0.3*X2 - 0.3*X3 + 0.6*X4 - 0.25*X5 + 0.3*X6 - 0.4*X1^2 + 0.3*X5^2
# f2 <- -0.3*X1 + 0.5*X2 + 0.2*X3 - 0.3*X4 - 0.5*X5 + 0.2*X6 + 0.3*X2^2 + 0.2*X6^2
# 
# f1 <- 0.5*X1 + 0.3*X2 - 0.3*X3 + 0.6*X4 - 0.25*X5 + 0.3*X6 - 0.2*X3*X4
# f2 <- -0.3*X1 + 0.5*X2 + 0.2*X3 - 0.3*X4 - 0.5*X5 + 0.2*X6 + 0.1*X3*X4
# 
# f1 <- 0.5*X1 + 0.3*X2 - 0.3*X3 + 0.6*X4 - 0.25*X5 + 0.3*X6 - 0.4*X1^2 - 0.2*X3*X4
# f2 <- -0.3*X1 + 0.5*X2 + 0.2*X3 - 0.3*X4 - 0.5*X5 + 0.2*X6 + 0.3*X2^2 + 0.1*X3*X4

simulate_multitreatment_data <- function(n = 1000, seed = 123) {
  set.seed(seed)
  
  # 1. Generate the four covariates:
  #    - X1, X2, X4 ~ N(0,1)
  #    - X3 ~ Bernoulli(0.5)
  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  X3 <- rnorm(n, mean = 0, sd = 1)
  X4 <- rnorm(n, mean = 0, sd = 1)
  
  # 2. Define linear predictors for the multinomial logistic model:
  #    T in {1, 2, 3}.
  #    Since X4 is *only* for outcome, we omit it from f1 and f2.
  
  f1 <-  0.5 * X1 + 0.3 * X2 + 0.4 * X3     # used for P(T=2)
  f2 <- -0.2 * X1 + 0.6 * X2 - 0.3 * X3     # used for P(T=3)
  
  # 3. Convert these into multinomial probabilities:
  #    p1 = 1 / [1 + exp(f1) + exp(f2)]
  #    p2 = exp(f1) / [1 + exp(f1) + exp(f2)]
  #    p3 = exp(f2) / [1 + exp(f1) + exp(f2)]
  
  denom <- 1 + exp(f1) + exp(f2)
  p1 <- 1         / denom
  p2 <- exp(f1)   / denom
  p3 <- exp(f2)   / denom
  
  # 4. Draw treatment T from {1, 2, 3} for each subject
  T <- numeric(n)
  for (i in seq_len(n)) {
    T[i] <- sample(x = 1:3, size = 1, prob = c(p1[i], p2[i], p3[i]))
  }
  
  # 5. Define the outcome model:
  #    E[Y|X] = 4 + 0.3*X1 - 0.2*X2 + 0.6*X4
  #            + (-0.5 * I{T=2}) + (0.7 * I{T=3})
  #    Then add Normal(0,1) noise.
  
  base_mean <- 4 + 0.3 * X1 - 0.2 * X2 + 0.6 * X4
  Y <- base_mean +
    ifelse(T == 2, -0.5, 0) +
    ifelse(T == 3,  0.7, 0) +
    rnorm(n, mean = 0, sd = 1)  # noise
  
  # Return a data frame
  return(data.frame(
    id        = 1:n,
    X1        = X1,
    X2        = X2,
    X3        = X3,
    X4        = X4,
    Treatment = T,
    Y         = Y
  ))
}
################################################################################################
# Estimate the generalized propensity scores
library(nnet)
estimate_propensity_scores <- function(data) {
  data$Treatment <- as.factor(data$Treatment)
  
  # Fit the multinomial logistic regression model
  model <- multinom(Treatment ~ . - Treatment - id - Y, data = data, trace = FALSE)
  
  # Predict probabilities for each treatment (propensity scores)
  propensity_scores <- predict(model, type = "probs")
  treatment_levels <- levels(data$Treatment)
  colnames(propensity_scores) <- paste0("PS", treatment_levels)
  return(propensity_scores)
}
################################################################################################
# Compute the common support set
common_support <- function(ps_matrix, treatment_col) {
  # Inputs:
  # - ps_matrix: Data frame or matrix where each row is a subject,
  #              columns are propensity scores for treatments,
  #              and one column specifies the treatment received.
  # - treatment_col: The name of the column indicating the treatment received.
  
  # Extract treatment groups
  treatment_groups <- unique(ps_matrix[[treatment_col]])
  
  # Initialize results
  low_values <- numeric(length(treatment_groups))
  high_values <- numeric(length(treatment_groups))
  
  # Loop through each treatment (target treatment)
  for (target_treatment in treatment_groups) {
    # Extract propensity scores for the target treatment (column corresponding to target_treatment)
    ps_target <- ps_matrix[[paste0("PS", target_treatment)]]
    
    # Compute min and max for each actual treatment group
    group_min <- sapply(treatment_groups, function(actual_treatment) {
      # Filter rows for subjects in the current actual treatment group
      group_data <- ps_matrix[ps_matrix[[treatment_col]] == actual_treatment, ]
      # Extract propensity scores for the target treatment
      ps_values <- group_data[[paste0("PS", target_treatment)]]
      return(min(ps_values))
    })
    
    group_max <- sapply(treatment_groups, function(actual_treatment) {
      # Filter rows for subjects in the current actual treatment group
      group_data <- ps_matrix[ps_matrix[[treatment_col]] == actual_treatment, ]
      # Extract propensity scores for the target treatment
      ps_values <- group_data[[paste0("PS", target_treatment)]]
      return(max(ps_values))
    })
    
    # Compute r(t, X)^(low) and r(t, X)^(high) for the current target treatment
    low_values[target_treatment] <- max(group_min)
    high_values[target_treatment] <- min(group_max)
  }
  
  # Return as a data frame for each treatment
  # The dimension for Low and High corresponds to the number of treatments.
  return(list(
    Low = low_values,
    High = high_values
  ))
}

# Function to filter subjects based on common support intervals
filter_valid_subjects_indices <- function(ps_data, css) {
  # It is checking whether each PS for each subject lies in the common support.
  # For example, c(2, 4, 6) is valid if low = c(1, 2, 3) and high = c(5, 6, 7)
  # as 2 in [1, 3], 4 in [2, 6], 6 in [3, 7]
  low <- css$Low
  high <- css$High
  
  # Check if each propensity score is within its corresponding [low, high] range
  is_valid <- apply(ps_data, 1, function(row) {
    all(row >= low & row <= high)
  })
  
  # Return the indices of valid rows
  valid_indices <- which(is_valid)
  return(valid_indices)
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

perform_matching <- function(data, ref_treatment, t_prime, epsilon = 0.25) {
  # Initialize a list to store matched data for each cluster
  all_matches <- data.frame()
  
  # Loop through each cluster
  unique_clusters <- unique(data$Cluster)
  for (k in unique_clusters) {
    # Subset data for the current cluster
    cluster_data <- data[data$Cluster == k, ]
    
    # Extract logit-transformed propensity scores for the treatments
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    
    # Prepare treatment indicator
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1, 
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    
    # Perform matching using Matching::Match
    # Create the caliper
    caliper_width <- epsilon * sd(cluster_data$logit_r_t, na.rm = TRUE)
    
    # Define variables for Matching
    X <- cluster_data$logit_r_t  # Covariate for matching
    Tr <- cluster_data$treatment_indicator  # Treatment indicator
    
    # Perform 1:1 nearest-neighbor matching with caliper
    match_result <- Match(
      Y = NULL,  # No outcome yet
      Tr = Tr,
      X = X,
      M = 1,  # 1:1 matching
      caliper = caliper_width,
      replace = TRUE  # With replacement
    )
    
    # Extract matched indices
    matched_id <- data.frame(id_ref = cluster_data$id[match_result$index.treated],
                             t_ref = cluster_data$Treatment[match_result$index.treated],
                             id_t_prime = cluster_data$id[match_result$index.control],
                             t_prime = cluster_data$Treatment[match_result$index.control])
    # cluster = k)
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################
# Matching for multiple matches such as 1:2
perform_matching_2 <- function(data, ref_treatment, t_prime, epsilon = 0.25, n_match) {
  # Initialize a list to store matched data for each cluster
  all_matches <- data.frame()
  
  # Loop through each cluster
  unique_clusters <- unique(data$Cluster)
  for (k in unique_clusters) {
    # Subset data for the current cluster
    cluster_data <- data[data$Cluster == k, ]
    
    # Extract logit-transformed propensity scores for the treatments
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    
    # Prepare treatment indicator
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1, 
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    
    # Perform matching using Matching::Match
    # Create the caliper
    caliper_width <- epsilon * sd(cluster_data$logit_r_t, na.rm = TRUE)
    
    # Define variables for Matching
    X <- cluster_data$logit_r_t  # Covariate for matching
    Tr <- cluster_data$treatment_indicator  # Treatment indicator
    
    # Perform 1:1 nearest-neighbor matching with caliper
    match_result <- Match(
      Y = NULL,  # No outcome yet
      Tr = Tr,
      X = X,
      M = n_match,  # 1:n_match matching
      caliper = caliper_width,
      replace = TRUE  # With replacement
    )
    
    # Extract matched indices
    matched_id <- data.frame(id_ref = cluster_data$id[match_result$index.treated],
                             t_ref = cluster_data$Treatment[match_result$index.treated],
                             id_t_prime = cluster_data$id[match_result$index.control],
                             t_prime = cluster_data$Treatment[match_result$index.control])
    # cluster = k)
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################
# VM_MD
perform_matching_vm_md <- function(data, ref_treatment, t_prime, epsilon = 0.25) {
  all_matches <- data.frame()
  
  unique_clusters <- unique(data$Cluster)
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    # We'll create two logit columns: logit_r_t & logit_r_tprime
    logit_transform <- function(p) log(p / (1 - p))
    
    cluster_data$logit_r_t      <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    cluster_data$logit_r_tprime <- logit_transform(cluster_data[[paste0("PS", t_prime)]])
    
    # Keep only those with T in {ref_treatment, t_prime}
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    # Construct the 2D matrix for each subject
    mat_2d <- as.matrix(cluster_data[, c("logit_r_t", "logit_r_tprime")])
    
    # Compute covariance & inverse-cholesky to get Mahalanobis
    Sigma <- cov(mat_2d)
    # If Sigma is singular or very small cluster, you might do tryCatch() or fallback
    Sigma_inv <- solve(Sigma)
    L <- chol(Sigma_inv)  # upper triangular
    
    # Transform each row: WeightedX_i = L %*% mat_2d_i
    WeightedX <- t(L %*% t(mat_2d))
    
    # For the caliper, we can define the caliper in WeightedX space
    WeightedX_mean <- apply(WeightedX, 1, function(z) sqrt(sum(z^2)))  # a typical distance magnitude
    caliper_width <- epsilon * sd(WeightedX_mean, na.rm = TRUE)
    
    # Now match on WeightedX (2D). Matching::Match will compute Euclidean distance in WeightedX
    # which is Mahalanobis distance in original logit space.
    match_result <- Match(
      Y = NULL,
      Tr = cluster_data$treatment_indicator,
      X = WeightedX,         # pass 2D
      M = 1,                 # 1-to-1 matching
      caliper = caliper_width,
      replace = TRUE
    )
    
    matched_id <- data.frame(
      id_ref     = cluster_data$id[match_result$index.treated],
      t_ref      = cluster_data$Treatment[match_result$index.treated],
      id_t_prime = cluster_data$id[match_result$index.control],
      t_prime    = cluster_data$Treatment[match_result$index.control],
      #cluster    = k
    )
    all_matches <- rbind(all_matches, matched_id)
  }
  return(all_matches)
}
################################################################################
# VM_MDnc
perform_matching_vm_mdnc <- function(data, ref_treatment, t_prime) {
  all_matches <- data.frame()
  
  unique_clusters <- unique(data$Cluster)
  for (k in unique_clusters) {
    cluster_data <- data[data$Cluster == k, ]
    
    # 1) Create 2D logit space
    logit_transform <- function(p) log(p / (1 - p))
    cluster_data$logit_r_t      <- logit_transform(cluster_data[[paste0("PS", ref_treatment)]])
    cluster_data$logit_r_tprime <- logit_transform(cluster_data[[paste0("PS", t_prime)]])
    
    # 2) Keep only T in {t, t'}
    cluster_data$treatment_indicator <- ifelse(
      cluster_data$Treatment == ref_treatment, 1,
      ifelse(cluster_data$Treatment == t_prime, 0, NA)
    )
    cluster_data <- subset(cluster_data, !is.na(treatment_indicator))
    if (nrow(cluster_data) < 2) next
    
    mat_2d <- as.matrix(cluster_data[, c("logit_r_t", "logit_r_tprime")])
    
    # 3) Cov & inverse-cholesky for Mahalanobis
    Sigma <- cov(mat_2d)
    Sigma_inv <- solve(Sigma)
    L <- chol(Sigma_inv)
    WeightedX <- t(L %*% t(mat_2d))
    
    # 4) "No caliper" => caliper = 0 
    match_result <- Match(
      Y = NULL,
      Tr = cluster_data$treatment_indicator,
      X = WeightedX,  # 2D
      M = 1,
      caliper = 0,    # No distance restriction
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
################################################################################
calc_standardized_bias <- function(data,
                                          covariates,
                                          ref_treatment = 1) {
  # 1) Number of matched sets (pairs):
  n_pairs <- length(unique(data$pair))
  
  # 2) psi_i: frequency of each subject's appearance
  count_by_id <- table(data$id)
  psi_lookup  <- setNames(as.numeric(count_by_id), names(count_by_id))
  
  get_psi <- function(id_i) psi_lookup[as.character(id_i)]
  
  # 3) Identify the two treatments actually present
  trt_levels <- sort(unique(data$Treatment))
  if (length(trt_levels) != 2) {
    stop("Data must have exactly 2 treatments!")
  }
  
  # 4) Weighted mean of each covariate by treatment
  wmean_mat <- matrix(NA,
                      nrow = length(covariates),
                      ncol = 2,
                      dimnames = list(covariates, paste0("T", trt_levels)))
  
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
      # Weighted mean = sum(...) / number_of_pairs
      wmean_mat[p, paste0("T", trt)] <- weighted_sum / n_pairs
    }
  }
  
  # 5) Weighted SD in the reference group
  wsd_vec <- numeric(length(covariates))
  names(wsd_vec) <- covariates
  
  ref_col <- paste0("T", ref_treatment)
  rows_ref <- which(data$Treatment == ref_treatment)
  if (length(rows_ref) == 0) {
    stop("No subjects found for reference treatment in the data!")
  }
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
  
  # 6) Standardized difference for the other treatment vs. ref_treatment
  other_treatment <- setdiff(trt_levels, ref_treatment)
  if (length(other_treatment) != 1) {
    stop("Could not identify exactly one 'other' treatment.")
  }
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
  
  # 7) Absolute SB and "Max2SB"
  #    With only 1 pair, max2SB_p = |SB_p|.
  abs_sb  <- abs(sb_vec)
  max2SBp <- abs_sb
  
  return(list(
    weighted_mean     = wmean_mat,  # Weighted means for each cov and treat
    sb_denominator    = wsd_vec,    # Weighted SD in reference group
    standardized_diff = sb_vec,     # SB for each covariate
    abs_sb            = abs_sb,     # Absolute SB
    max2SB            = max2SBp     # The "Max2SB" for each covariate
  ))
}


rrrrr <- calc_standardized_bias(
  data         = final_cohort_1_2,
  covariates   = c("X1","X2","X3","X4"),
  ref_treatment = 1
)

################################################################################
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

################################################################################
# Example usage
data <- simulate_multitreatment_data(3000)

# Step 1: Estimate the generalized propensity scores
ps <- estimate_propensity_scores(data)
ps_trt <- as.data.frame(cbind(ps, Treatment = data$Treatment))

# Step 2: Drop the invalid subjects
## Calculate the common support set
css <- common_support(ps_trt, "Treatment")
## Find the valid subjects
valid_index <- filter_valid_subjects_indices(ps, css)
valid_data <- data[valid_index,]
## Drop and compute GPS again
ps_new <- estimate_propensity_scores(valid_data)
valid_data_ps <- cbind(valid_data, ps_new)

# Step 3: KMC and matching
## For the first pair
### KMC
kmc_A <- perform_kmc(ps_new, ref_treatment = "PS1", t_prime = "PS2", K = 3)
valid_data_ps_A <- cbind(valid_data_ps, Cluster = kmc_A$clusters)
### Matching within each cluster
# matches_1_2 <- perform_matching(valid_data_ps_A, ref_treatment = 1, t_prime = 2)
# matches_1_2 <- perform_matching_2(valid_data_ps_A, ref_treatment = 1, t_prime = 2, n_match = 2)
matches_1_2 <- perform_matching_vm_md(valid_data_ps_A, ref_treatment = 1, t_prime = 2)
## For the second pair
### KMC
kmc_B <- perform_kmc(ps_new, ref_treatment = "PS1", t_prime = "PS3", K = 3)
valid_data_ps_B <- cbind(valid_data_ps, Cluster = kmc_B$clusters)
### Matching within each cluster
matches_1_3 <- perform_matching(valid_data_ps_B, ref_treatment = 1, t_prime = 3)
# matches_1_3 <- perform_matching_2(valid_data_ps_B, ref_treatment = 1, t_prime = 3, nmatch = 2)

kmc_C <- perform_kmc(ps_new, ref_treatment = "PS2", t_prime = "PS3", K = 3)
valid_data_ps_C <- cbind(valid_data_ps, Cluster = kmc_C$clusters)
### Matching within each cluster
matches_2_3 <- perform_matching(valid_data_ps_C, ref_treatment = 2, t_prime = 3)
# matches_2_3 <- perform_matching_2(valid_data_ps_C, ref_treatment = 2, t_prime = 3, nmatch = 2)

# Step 4: Composite to the final matched cohort.
final_cohort_1_2 <- data.frame(
  id   = c(rbind(matches_1_2$id_ref, matches_1_2$id_t_prime)),
  t    = c(rbind(matches_1_2$t_ref, matches_1_2$t_prime)),
  pair = rep(seq_len(nrow(matches_1_2)), each = 2)
)

final_cohort_1_2 <- left_join(
  final_cohort_1_2,
  valid_data_ps, 
  by = c("id", "t" = "Treatment")
) %>%
  rename(Treatment = t)

# === 2) Construct the final matched cohort for the (1 vs 3) pairs ===
final_cohort_1_3 <- data.frame(
  id   = c(rbind(matches_1_3$id_ref, matches_1_3$id_t_prime)),
  Treatment    = c(rbind(matches_1_3$t_ref, matches_1_3$t_prime)),
  pair = rep(seq_len(nrow(matches_1_3)), each = 2)
)

final_cohort_1_3 <- left_join(
  final_cohort_1_3,
  valid_data_ps,
  by = c("id", "Treatment")
)

# === 3) Construct the final matched cohort for the (2 vs 3) pairs ===
final_cohort_2_3 <- data.frame(
  id   = c(rbind(matches_2_3$id_ref, matches_2_3$id_t_prime)),
  Treatment    = c(rbind(matches_2_3$t_ref, matches_2_3$t_prime)),
  pair = rep(seq_len(nrow(matches_2_3)), each = 2)
)

final_cohort_2_3 <- left_join(
  final_cohort_2_3,
  valid_data_ps,
  by = c("id", "Treatment")
)
################################################################################