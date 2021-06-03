log_odds <- function(p) log(p[1]) - log(p[2]) - log(p[3]) + log(p[4])

log_omega_estimator <- function(counts, lambda)
{
  count_00_smoothed <- (1 - lambda)^2 * counts[1] + lambda * (1 - lambda) * (counts[2] + counts[3]) + lambda^2 * counts[4]
  count_01_smoothed <- (1 - lambda)^2 * counts[2] + lambda * (1 - lambda) * (counts[1] + counts[4]) + lambda^2 * counts[3]
  count_10_smoothed <- (1 - lambda)^2 * counts[3] + lambda * (1 - lambda) * (counts[1] + counts[4]) + lambda^2 * counts[2]
  count_11_smoothed <- (1 - lambda)^2 * counts[4] + lambda * (1 - lambda) * (counts[2] + counts[3]) + lambda^2 * counts[1]
  counts_smoothed <- c(count_00_smoothed, count_01_smoothed, count_10_smoothed, count_11_smoothed)
  log_omega_ <- log_odds(counts_smoothed)
  
  return(log_omega_)
}

optimal_lambda <- function(n, p)
{
  a <- (p[2] + p[3]) * (1/p[1] + 1/p[4]) - (p[1] + p[4]) * (1/p[2] + 1/p[3])
  b <- 0.5 * (1/p[1] - 1/p[2] - 1/p[3] + 1/p[4])
  k00 <- 1/p[2] + 1/p[3] + 2/p[1]
  k01 <- 1/p[1] + 1/p[4] + 2/p[2]
  k10 <- 1/p[1] + 1/p[4] + 2/p[3]
  k11 <- 1/p[2] + 1/p[3] + 2/p[4]
  
  num <- a * b + k00 + k01 + k10 + k11
  den <- n * a^2 + k00^2 * p[1] + k01^2 * p[2] + k10^2 * p[3] + k11^2 * p[4]
  optimal_lambda_ <- num / den
  
  return(optimal_lambda_)
}

optimal_lambda_pi <- function(counts, c=0.5)
{
  n_pi <- sum(counts) + 4 * c
  p_pi <- (counts + c) / n_pi
  optimal_lambda_pi_ <- optimal_lambda(n_pi, p_pi)
  return(optimal_lambda_pi_)
}

log_omega_estimator_pi <- function(counts, c=0.5)
{
  optima_lambda_pi_ <- optimal_lambda_pi(counts, c)
  log_omega_ <- log_omega_estimator(counts, optima_lambda_pi_)
  return(log_omega_)
}

log_omega_estimator_simple <- function(counts, c=0)
{
  counts_smoothed <- counts + c
  log_omega_ <- log_odds(counts_smoothed)
  return(log_omega_)
}