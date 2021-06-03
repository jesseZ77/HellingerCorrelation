library(tidyverse)
source("util.R")
source("mse.R")
source("estimator.R")
# setwd("Simulations/")

m <- 5000
n_list <- c(seq(50, 200, 5), seq(250, 500, 10), seq(550, 1000, 50), seq(1000, 2000, 100))

omega_list <- c(2, 10, 20, 50)
p0x_list <- c(0.1, 0.33, 0.66, 0.9)
px0_list <- c(0.1, 0.33, 0.66, 0.9)

p_mat_list <- expand.grid(omega_list, p0x_list, px0_list)

ctr <- 1
max_ctr <- dim(p_mat_list)[1] * length(n_list)
output <- matrix(0L, max_ctr, 9)

for(i in 1:dim(p_mat_list)[1])
{
  omega <- p_mat_list[i, 1]
  p0x <- p_mat_list[i, 2]
  px0 <- p_mat_list[i, 3]
  p <- generate_p_matrix(omega, p0x, px0)

  for(n in n_list)
  {
    lambda_opt <- optimal_lambda(n, p)
    sample_ <- generate_m_samples(m, n, p)
    estimates_ <- apply(sample_, 2, function(x) log_omega_estimator(x, lambda_opt))

    bias_e <- mean(estimates_) - log(omega)
    var_e <- var(estimates_)
    bias_a <- bias(n, p, lambda_opt)
    var_a <- variance(n, p, lambda_opt)

    output[ctr, ] <- c(omega, p0x, px0, n, lambda_opt, bias_e, var_e, bias_a, var_a)

    if(ctr %% 1000 == 0) print(paste("Done", ctr, "out of", max_ctr))
    ctr <- ctr + 1
  }
}
print("Finished")

df_results <- data.frame(output)
colnames(df_results) <- c("omega", "p0x", "px0", "n", "lambda", "bias_e", "var_e", "bias_a", "var_a")

df_results$bias_sq_e <- df_results$bias_e^2
df_results$bias_sq_a <- df_results$bias_a^2

df_results$mse_e <- df_results$bias_sq_e + df_results$var_e
df_results$mse_a <- df_results$bias_sq_a + df_results$var_a

df_results$mse_diff <- df_results$mse_a - df_results$mse_e
df_results$bias_sq_diff <- df_results$bias_sq_a - df_results$bias_sq_e
df_results$var_diff <- df_results$var_a - df_results$var_e

write.csv(df_results, file="results_mse.csv", row.names = FALSE)

df_results <- read.csv("results_mse.csv")

for (omega_fixed in omega_list)
{
  df_plot <- df_results %>% 
    filter(omega == omega_fixed, n >= 250, n <= 2000)
  
  p_mse <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=mse_e, colour="Empirical")) +
    geom_line(aes(y=mse_a, colour="Analytical")) +
    scale_colour_manual("Legend", values = c("blue", "green")) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("MSE for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p_bias_sq <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=bias_sq_e, colour="Empirical")) +
    geom_line(aes(y=bias_sq_a, colour="Analytical")) +
    scale_colour_manual("Legend", values = c("blue", "green")) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("Bias squared for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p_var <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=var_e, colour="Empirical")) +
    geom_line(aes(y=var_a, colour="Analytical")) +
    scale_colour_manual("Legend", values = c("blue", "green")) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("Variance for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0("plots/mse_plot_", omega_fixed, ".png"), p_mse, width=12, height=12)
  ggsave(paste0("plots/bias_sq_plot_", omega_fixed, ".png"), p_bias_sq, width=12, height=12)
  ggsave(paste0("plots/var_plot_", omega_fixed, ".png"), p_var, width=12, height=12)
  
  df_plot <- df_results %>% 
    filter(omega == omega_fixed, n >= 100, n <= 2000)
  
  p_mse_diff <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=mse_diff)) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("MSE diff for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p_bias_sq_diff <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=bias_sq_diff)) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("Bias squared diff for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p_var_diff <- ggplot(df_plot, aes(x=n)) +
    geom_line(aes(y=var_diff)) +
    facet_wrap(~p0x + px0, nrow=4, scales="free") + 
    ggtitle(paste("Variance diff for omega = ", omega_fixed)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0("plots/mse_diff_plot_", omega_fixed, ".png"), p_mse_diff, width=12, height=12)
  ggsave(paste0("plots/bias_sq_diff_plot_", omega_fixed, ".png"), p_bias_sq_diff, width=12, height=12)
  ggsave(paste0("plots/var_diff_plot_", omega_fixed, ".png"), p_var_diff, width=12, height=12)
}

