generate_p_matrix <- function(omega, p0x, px0)
{
  a <- (omega - 1)
  b <- (1 - omega) * p0x + (1 - omega) * px0 - 1
  c <- omega * p0x * px0
  p00 <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  p01 <- p0x - p00
  p10 <- px0 - p00
  p11 <- 1 - p00 - p01 - p10
  p <- c(p00, p01, p10, p11)
  
  return(p)
}

generate_m_samples <- function(m, n, p) replicate(m, rpois(4, n * p))