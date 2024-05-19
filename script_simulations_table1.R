# Libraries ----
library(tidyverse)

# Simulations ----

## Parameters ----

### Means ----
mu_m_1 <- 1
mu_m_2 <- 2

### Variance
sigma2_m_1 <- 1

### Null hypothesis ----
null_hypothesis <- mu_m_1 - mu_m_2

### Type I error rates
alpha <- c(0.1, 0.05)

### Replications ----
R <- 10000

### Proportions between the variances ----
rho <- seq.int(from = 2, to = 10, by = 2) |> 
  as.integer()

### Size of the 2, m_1 and m_2, samples ----
n <- 2:12

### Degrees of freedom Welch–Satterthwaite approach ---- 
df_ws <- function(n, rho) {
  
  df <- ((n - 1)*(rho + 1)^2) / (rho^2 + 1)
  
  return(df)
}

### t Welch–Satterthwaite critical values
t_cv_lower <- function(alpha, df) {
  
  lower <- qt(p = alpha/2, 
              df = df,
              lower.tail = TRUE)
  
  return(lower)
}

### Probability density function m_1 and m_2
probdf <- function(n, mean, sd) {
  
  return(rnorm(n, mean = mean, sd = sd))
}

### t-statistic
t_stat <- function(mu_1, mu_2, 
                   sigma2_1, sigma2_2,
                   n_1, n_2,
                   null_hypothesis) {
  
  statistic <- ((mu_1 - mu_2) - null_hypothesis) / sqrt(x = (sigma2_1/n_1) + (sigma2_2/n_2))
  
  return(statistic)
  
}

# Lower bound
d <- 0.10

# Fixing seed
set.seed(seed = 1234, 
         kind = "Mersenne-Twister", 
         normal.kind = "Inversion", 
         sample.kind = "Rejection")

simulation_tbl <- expand_grid(rho = rho, 
                              n = n,
                              R = 1:R) |> 
  mutate(samples_m_1 = map2(.x = R, 
                            .y = n,
                            .f = ~ probdf(n = .y,
                                          mean = mu_m_1,
                                          sd = sqrt(x = sigma2_m_1))),
         samples_m_2 = pmap(.l = list(R, n, rho),
                            .f = ~ probdf(n = ..2,
                                          mean = mu_m_2,
                                          sd = sqrt(x = ..3*sigma2_m_1)))) |> 
  mutate(samples_mu_m_1 = map(.x = samples_m_1, 
                              .f = mean),
         samples_mu_m_2 = map(.x = samples_m_2, 
                              .f = mean),
         samples_sigma2_m_1 = map(.x = samples_m_1, 
                                  .f = var),
         samples_sigma2_m_2 = map(.x = samples_m_2, 
                                  .f = var),
         .keep = "unused") |> 
  unnest(cols = c(samples_mu_m_1:samples_sigma2_m_2)) |> 
  mutate(t = t_stat(mu_1 = samples_mu_m_1,
                    mu_2 = samples_mu_m_2, 
                    sigma2_1 = samples_sigma2_m_1, 
                    sigma2_2 = samples_sigma2_m_2, 
                    n_1 = n, n_2 = n, 
                    null_hypothesis = null_hypothesis),
         rho_est = samples_sigma2_m_2/samples_sigma2_m_1) |> 
  select(-c(samples_mu_m_1:samples_sigma2_m_2)) 

simulation_tbl_stacked <- rep.int(x = list(simulation_tbl), 
                                  times = length(alpha)) |> 
  bind_rows() |> 
  mutate(alpha = rep(x = alpha, 
                     each = nrow(simulation_tbl))) |>
  mutate(df_ws = df_ws(rho = rho_est,
                       n = n),
         t_ws_cv_lower = t_cv_lower(alpha = alpha,
                                    df = df_ws)) |> 
  group_by(rho, n, alpha) |> 
  mutate(alpha_est = ecdf(x = t)(t_ws_cv_lower) + (1 - ecdf(x = t)(-t_ws_cv_lower))) |> 
  ungroup() |> 
  mutate(pct_rel_dist = (abs(alpha - alpha_est) / alpha) > d) |> 
  group_by( rho, n, alpha) |> 
  summarize(pct_greater_d = sum(pct_rel_dist) / n(),
            .groups = "drop")

# Export data ----
simulation_tbl_stacked |> 
  write_csv(file = "data/simulation_tbl_stacked.csv")
