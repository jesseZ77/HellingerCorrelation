bias <- function(n, p, lambda)
{
  a <- (p[2] + p[3]) * (1/p[1] + 1/p[4]) - (p[1] + p[4]) * (1/p[2] + 1/p[3])
  b <- 0.5 * (1/p[1] - 1/p[2] - 1/p[3] + 1/p[4])
  bias_ <- a * lambda - b / n
  
  return(bias_)
}

variance <- function(n, p, lambda)
{
  k00 <- 1/p[2] + 1/p[3] + 2/p[1]
  k01 <- 1/p[1] + 1/p[4] + 2/p[2]
  k10 <- 1/p[1] + 1/p[4] + 2/p[3]
  k11 <- 1/p[2] + 1/p[3] + 2/p[4]
  
  term00 <- p[1]/n * (1/p[1] - lambda * k00)^2
  term01 <- p[2]/n * (1/p[2] - lambda * k01)^2
  term10 <- p[3]/n * (1/p[3] - lambda * k10)^2
  term11 <- p[4]/n * (1/p[4] - lambda * k11)^2
  variance_ <- term00 + term01 + term10 + term11
  
  return(variance_)
}

mse <- function(n, p, lambda)
{
  bias_ <- bias(n, p, lambda)
  variance_ <- variance(n, p, lambda)
  mse_ <- bias_^2 + variance_
  
  return(mse_)
}
