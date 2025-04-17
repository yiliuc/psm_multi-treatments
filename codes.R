simulate_multitreatment_data <- function(nt1, gamma, tau, sigma2, sigma3, b, p) {
  library(MASS)  # for mvrnorm
  
  # Sample sizes
  nt2 <- gamma * nt1
  nt3 <- gamma^2 * nt1
  N <- nt1 + nt2 + nt3
  
  # Check p is divisible by 3
  if (p %% 3 != 0) stop("p must be divisible by 3")
  
  # Function to construct mean vector
  construct_mu <- function(pos, b, p) {
    if (p %% 3 != 0) stop("p must be a multiple of 3")
    mu <- numeric(p)
    repeats <- p / 3
    for (i in 0:(repeats - 1)) {
      mu[3 * i + pos] <- b
    }
    return(mu)
  }
  
  ########################################
  # Mean vectors
  mu1 <- construct_mu(1, b, p)
  mu2 <- construct_mu(2, 2, p)
  mu3 <- construct_mu(3, 3, p)
  ########################################
  
  # Covariance matrices
  make_sigma <- function(diag_val, tau, p) {
    mat <- matrix(tau, nrow = p, ncol = p)
    diag(mat) <- diag_val
    return(mat)
  }
  
  Sigma1 <- make_sigma(1, tau, p)
  Sigma2 <- make_sigma(sigma2, tau, p)
  Sigma3 <- make_sigma(sigma3, tau, p)
  
  # Generate X from multivariate normal
  X1 <- mvrnorm(n = nt1, mu = mu1, Sigma = Sigma1)
  X2 <- mvrnorm(n = nt2, mu = mu2, Sigma = Sigma2)
  X3 <- mvrnorm(n = nt3, mu = mu3, Sigma = Sigma3)
  
  # Combine all
  X <- rbind(X1, X2, X3)
  Tvec <- c(rep(1, nt1), rep(2, nt2), rep(3, nt3))
  
  # Return as a data.frame
  df <- data.frame(Treatment = Tvec, X)
  return(df)
}
################################################################################
# Generate the outcome, assume there is no intercept.
generate_outcomes <- function(df, betas, sigma_eps = 1) {
  if (!"Treatment" %in% names(df)) stop("Data frame must contain 'Treatment' column.")
  
  X <- as.matrix(df[, -1])  # Exclude Treatment column
  n <- nrow(X)
  p <- ncol(X)
  
  # Sanity check
  # if (length(alphas) != 3 || length(betas) != 3) {
  #   stop("alphas and betas must be lists of length 3")
  # }
  
  # Allocate Y
  Y <- numeric(n)
  
  for (i in 1:n) {
    t <- df$Treatment[i]
    # alpha_t <- alphas[[t]]
    beta_t <- betas[[t]]
    Y[i] <- X[i, ] %*% beta_t + rnorm(1, mean = 0, sd = sigma_eps)
  }
  
  df$Y <- Y
  return(df)
}
################################################################################
estimate_ATT_t1_t2 <- function(df, t1 = 1, t2 = 2) {
  # Compute psi_i: number of times unit appears in matched sets
  df$psi <- ave(df$pair, df$id, FUN = length)
  
  # n_trip: number of unique matched triplets
  n_trip <- length(unique(df$pair))
  
  # Indicator variables
  df$I_t1 <- as.numeric(df$Treatment == t1)
  df$I_t2 <- as.numeric(df$Treatment == t2)
  
  # Contribution to SATT
  df$weighted_contribution <- df$psi * (df$Y * df$I_t1 - df$Y * df$I_t2)
  
  # Estimate
  SATT_est <- sum(df$weighted_contribution) / n_trip
  return(SATT_est)
}

# data <- simulate_multitreatment_data(1000, 1, 0, 2, 3, 1, 3)
# betas <- list(rep(0.5, 3), rep(0.5, 3), rep(0.5, 3))
# data <- generate_outcomes(data, betas)
# data$id <- 1:nrow(data)
# 
# estimate_ATT_t1_t2(final_cohort, 1, 2)


# WE SHOULD INCLUDE THE TREATMENT IN GENERATING THE OUTCOME!!!!!!!!!!



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
