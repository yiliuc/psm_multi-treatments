################################################################################################
# Generate the covariates
generate_data <- function(sample_sizes, mu, sigma, tau, p) {
  # sample_sizes: the vector containing the sample sizes for each treatment
  # mu: the vector containing the mean value of covariates for each treatment
  # sigma: the vector containing the variance of covariates for each treatment
  # tau: the covariance
  # p: the number of variables
  
  # Note that I assume the treatment i has mean vector (0, ..., mu_i, ...., 0) and 
  # covariance matrix with diagonal (sigma_i, ..., sigma_i). That is why the 
  # input mu and sigma are both vectors. I assume the mean and variance for all
  # covariates in the same treatment are same!
  
  library(MASS)  # For multivariate normal generation
  
  T <- length(sample_sizes)  # Number of treatment groups
  # Generate data for each treatment group
  data_list <- list()
  
  for (t in 1:T) {
    # Define mean vector for treatment t
    mu_t <- rep(0, p)
    if (t <= p) {
      mu_t[t] <- mu[t]  # Assign mean only to the corresponding dimension
    }
    
    # Define covariance matrix for treatment t
    Sigma_t <- matrix(tau, nrow = p, ncol = p)
    diag(Sigma_t) <- sigma[t]  # Set the diagonal variances
    
    # Generate data for treatment t
    X_t <- mvrnorm(n = sample_sizes[t], mu = mu_t, Sigma = Sigma_t)
    
    # Store data with treatment label
    data_list[[t]] <- cbind(X_t, Treatment = t)
  }
  
  # Combine all treatment groups into a single dataframe
  data <- do.call(rbind, data_list)
  data <- as.data.frame(data)
  colnames(data)[1:p] <- paste0("X", 1:p)
  data$id <- 1:nrow(data)
  return(data)
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
################################################################################################
plot_ps_histogram <- function(data, ps, treatment1, treatment2, xlim, ylim, title) {
  t1_ps <- data[[ps]][data$Treatment == treatment1]
  t2_ps <- data[[ps]][data$Treatment == treatment2]
  
  # Create histograms
  hist(t1_ps, 
       breaks = 20,
       col = rgb(1, 0, 0, 0.5),
       xlim = xlim,
       ylim = ylim,
       ylab = "Count",
       main = title)
  hist(t2_ps, 
       breaks = 20,
       col = rgb(0, 0, 1, 0.5),
       add = TRUE)
}

plot_ps_histogram <- function(data, ps, treatment1, treatment2, xlim, ylim, title) {
  # A named vector: the names (i.e. "1","2","3") map to an RGBA color.
  # Adjust the colors as you like:
  color_map <- c(
    "1" = rgb(1, 0, 0, 0.5),  # semi-transparent red
    "2" = rgb(0, 1, 0, 0.5),  # semi-transparent green
    "3" = rgb(0, 0, 1, 0.5)   # semi-transparent blue
  )
  
  t1_ps <- data[[ps]][data$Treatment == treatment1]
  t2_ps <- data[[ps]][data$Treatment == treatment2]
  
  # Histogram for the first treatment group
  hist(t1_ps,
       breaks = 20,
       col    = color_map[as.character(treatment1)],
       xlim   = xlim,
       ylim   = ylim,
       ylab   = "Count",
       main   = title)
  
  # Histogram for the second treatment group, overlaid
  hist(t2_ps,
       breaks = 20,
       col    = color_map[as.character(treatment2)],
       add    = TRUE)
}

