library(tidyverse)
source("util.R")
source("mse.R")
source("estimator.R")
options(pillar.sigfig = 6)
# setwd("Simulations/")

m <- 5000
n_list <- c(500, 1000, 2000, 5000, 10000)

omega_list <- c(2, 10, 20, 50)
p0x_list <- c(0.1, 0.5, 0.66)
px0_list <- c(0.33, 0.5, 0.9)

# omega_list <- c(2)
# p0x_list <- c(0.1)
# px0_list <- c(0.1)

p_mat_list <- expand.grid(omega_list, p0x_list, px0_list)

# df_results <- tibble(
#   omega=numeric(0), p0x=numeric(0), px0=numeric(0),
#   lambda_opt=numeric(0), lambda_opt_pi=numeric(0),
#   estimate_oracle=numeric(0), estimate_pi=numeric(0), 
#   estimate_mle=numeric(0), estimate_c05=numeric(0)
# )

df_results <- tibble()
ctr <- 1
max_ctr <- dim(p_mat_list)[1] * length(n_list)

for(i in 1:dim(p_mat_list)[1])
{
  omega <- p_mat_list[i, 1]
  p0x <- p_mat_list[i, 2]
  px0 <- p_mat_list[i, 3]
  p <- generate_p_matrix(omega, p0x, px0)

  for(n in n_list)
  {
    sample_ <- generate_m_samples(m, n, p)
    
    lambda_opt <- optimal_lambda(n, p)
    estimate_oracle <- apply(sample_, 2, log_omega_estimator, lambda_opt)
    
    lambda_opt_pi <- apply(sample_, 2, optimal_lambda_pi, 0.5)
    estimate_pi <- apply(sample_, 2, log_omega_estimator_pi, 0.5)
    
    estimate_mle <- apply(sample_, 2, log_omega_estimator_c, 0)
    estimate_c05 <- apply(sample_, 2, log_omega_estimator_c, 0.5)
    
    df_tmp <- tibble(omega, p0x, px0, n,
                     lambda_opt, lambda_opt_pi,
                     estimate_oracle, estimate_pi, estimate_mle, estimate_c05)
    
    df_results <- rbind(df_results, df_tmp)
    ctr <- ctr + 1
    if(ctr %% 10 == 0) print(paste("Done", ctr, "out of", max_ctr))
  }
}
print("Finished")

write.csv(df_results, file="results_comparison.csv", row.names = FALSE)

df_plot <- df_results %>% 
  select(-c(lambda_opt, lambda_opt_pi)) %>% 
  pivot_longer(cols=c(estimate_oracle, estimate_pi, estimate_mle, estimate_c05), 
               names_to="estimator", values_to="estimate")
df_plot$estimate[df_plot$estimate == Inf] <- NaN
omega_fixed <- 10
p0x_fixed <- 0.1
px0_fixed <- 0.9

df_plot_tmp <- df_plot %>% filter(omega == omega_fixed, p0x == p0x_fixed, px0 == px0_fixed)

p_boxplot <- ggplot(df_plot_tmp, aes(x=estimator, y=estimate)) + 
  geom_boxplot() + 
  geom_hline(yintercept=log(omega_fixed), color="red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(~n + p0x + px0, nrow=1)

ggsave("plots/comparison_boxplot.png", p_boxplot, width=12, height=12)

df_summary <- df_plot %>% 
  group_by(omega, p0x, px0, n, estimator) %>% 
  summarise(mean=mean(estimate, na.rm=T), var=var(estimate, na.rm=T)) %>% 
  mutate(bias=mean - log(omega), bias_sq=bias^2, mse=bias_sq + var)

write.csv(df_summary, file="results_comparison_summary.csv", row.names = FALSE)

omega_fixed <- 50
p0x_fixed <- 0.1
px0_fixed <- 0.9

# df_summary_mse <- df_summary %>%
#   pivot_wider(id_cols=c(omega, p0x, px0, n), names_from=estimator, values_from=mse)

df_plot2 <- df_summary %>% filter(p0x == p0x_fixed, px0 == px0_fixed)

df_plot2 %>% 
  pivot_wider(id_cols=c(omega, p0x, px0, n), names_from=estimator, values_from=mse) %>% 
  filter(omega == 10)

ggplot(df_plot2, aes(x=estimator, y=mse)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~omega + n, scales="free_y")


