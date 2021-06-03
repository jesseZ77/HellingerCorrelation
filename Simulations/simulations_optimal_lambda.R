library(tidyverse)
source("util.R")
source("mse.R")
source("estimator.R")
# setwd("Simulations/")

m <- 5000
n_list <- c(seq(50, 500, 50), seq(600, 1500, 100), seq(2000, 5000, 500), seq(6000, 10000, 1000), seq(12000, 20000, 2000))

omega <- 10
p0x <- 0.3
px0 <- 0.6
p <- generate_p_matrix(omega, p0x, px0)

optimal_lambda_e <- rep(0, length(n_list))
optimal_lambda_a <- rep(0, length(n_list))
df_mse <- tibble(n=numeric(0), lambda=numeric(0), mse=numeric(0))

for(i in 1:length(n_list))
{
  n <- n_list[i]

  sample_ <- generate_m_samples(m, n, p)
  lambda_opt_a <- optimal_lambda(n, p)
  optimal_lambda_a[i] <- lambda_opt_a

  lambda_list <- seq(0, min(5 * lambda_opt_a, 0.2), length.out=100)
  mse_list <- rep(0, length(lambda_list))

  for(j in 1:length(lambda_list))
  {
    lambda <- lambda_list[j]
    estimates_ <- apply(sample_, 2, function(x) log_omega_estimator(x, lambda))
    bias_e <- mean(estimates_) - log(omega)
    var_e <- var(estimates_)
    mse_e <- bias_e^2 + var_e
    mse_list[j] <- mse_e
  }

  lambda_opt_e <- lambda_list[which.min(mse_list)]
  optimal_lambda_e[i] <- lambda_opt_e

  df_mse_tmp <- tibble(n=n, lambda=lambda_list, mse=mse_list)
  df_mse <- rbind(df_mse, df_mse_tmp)

  print(paste("Done", i, "out of", length(n_list)))
}

df_results <- tibble(n=n_list, lambda_a=optimal_lambda_a, lambda_e=optimal_lambda_e)
df_results$lambda_diff <- df_results$lambda_a - df_results$lambda_e

write.csv(df_results, file="results_optimal_lambda.csv", row.names = FALSE)
write.csv(df_mse, file="results_optimal_lambda_mse.csv", row.names = FALSE)

df_results <- read.csv("results_optimal_lambda.csv")

df_plot <- df_results %>% filter(n >= 250, n <= 20000)
p_comparison <- ggplot(df_plot, aes(x=n)) + 
  geom_line(aes(y=lambda_e, colour="Empirical")) +
  geom_line(aes(y=lambda_a, colour="Theoretical")) +
  scale_colour_manual("Legend", values = c("blue", "green")) +
  ggtitle("Optimal lambda vs n") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("n") +
  ylab("lambda")

p_diff <- ggplot(df_plot, aes(x=n, y=lambda_diff)) + 
  geom_line() + 
  ggtitle("Lambda diff vs n (analytical - empirical)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("n") +
  ylab("Diff")

ggsave("plots/optimal_lambda_comparison.png", p_comparison, width=8, height=8)
ggsave("plots/optimal_lambda_diff.png", p_diff, width=8, height=8)
