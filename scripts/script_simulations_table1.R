# Libraries ----
library(tidyverse)

# Function ----

ws_simulation <- function(mu_m_1 = 1,
                          mu_m_2 = 2,
                          sigma2_m_1 = 1,
                          rho = as.integer(seq.int(from = 2, to = 12, by = 2)),
                          n_1 = 2:14,
                          k = 1,
                          alpha = c(0.1, 0.05),
                          R = 10000,
                          d = 0.10,
                          seed = 1234) {
  
  # Parameters ----
  
  ## Means ----
  mu_m_1 <- mu_m_1
  mu_m_2 <- mu_m_2
  
  ## Variance ----
  sigma2_m_1 <- sigma2_m_1
  
  ### Proportions between the variance: sigma2_m_2 / sigma2_m_1 = rho
  rho <- rho 
  
  ## Sample size ----
  
  ### Size of the m_1 samples
  n_1 <- n_1
  
  ### Proportion between the samples size: n_2 / n_1 = k
  k <- k
  
  ## Degrees of freedom Welch–Satterthwaite approach ----
  ### nu in the paper
  df_ws <- function(n_1, rho, k) {
    
    df <- ((n_1 - 1)*(k*n_1 - 1)*(k + rho)^2) / ((k*n_1 - 1)*k^2 + (n_1 - 1)*rho^2)
    
    return(df)
  }
  
  ## Statistic ----
  
  ### Null hypothesis
  null_hypothesis <- mu_m_1 - mu_m_2
  
  ### t-statistic
  t_stat <- function(mu_1, mu_2, 
                     sigma2_1, sigma2_2,
                     n_1, n_2,
                     null_hypothesis) {
    
    statistic <- ((mu_1 - mu_2) - null_hypothesis) / sqrt(x = (sigma2_1/n_1) + (sigma2_2/n_2))
    
    return(statistic)
    
  }
  
  ## t Welch–Satterthwaite ----
  
  ### Type I error rates
  alpha <- alpha
  
  ### Critical values
  t_cv_lower <- function(alpha, df) {
    
    lower <- qt(p = alpha/2, 
                df = df,
                lower.tail = TRUE)
    
    return(lower)
  }
  
  ## Replications ----
  R <- R
  
  ## Probability density function m_1 and m_2 ----
  probdf <- function(n, mean, sd) {
    
    return(rnorm(n, mean = mean, sd = sd))
  }
  
  ## Lower bound ----
  d <- d
  
  ## Seed ----
  set.seed(seed = seed, 
           kind = "Mersenne-Twister", 
           normal.kind = "Inversion", 
           sample.kind = "Rejection")
  
  # Simulations ----
  simulation_tbl <- expand_grid(rho = rho, 
                                n_1 = n_1,
                                R = 1:R) |> 
    mutate(samples_m_1 = map2(.x = R, 
                              .y = n_1,
                              .f = ~ probdf(n = .y,
                                            mean = mu_m_1,
                                            sd = sqrt(x = sigma2_m_1))),
           samples_m_2 = pmap(.l = list(R, n_1, rho),
                              .f = ~ probdf(n = ..2*k,
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
                      n_1 = n_1, n_2 = k*n_1, 
                      null_hypothesis = null_hypothesis),
           rho_est = samples_sigma2_m_2/samples_sigma2_m_1) |> 
    select(-c(samples_mu_m_1:samples_sigma2_m_2))
  
  simulation_tbl_stacked <- rep.int(x = list(simulation_tbl), 
                                    times = length(alpha)) |> 
    bind_rows() |> 
    mutate(alpha = rep(x = alpha, 
                       each = nrow(simulation_tbl))) |>
    mutate(df_ws = df_ws(n_1 = n_1,
                         rho = rho_est,
                         k = k),
           t_ws_cv_lower = t_cv_lower(alpha = alpha,
                                      df = df_ws)) |> 
    group_by(rho, n_1, alpha) |> 
    mutate(alpha_est = ecdf(x = t)(t_ws_cv_lower) + (1 - ecdf(x = t)(-t_ws_cv_lower))) |> 
    ungroup() |> 
    mutate(pct_rel_dist = (abs(alpha - alpha_est) / alpha) > d) |> 
    group_by( rho, n_1, alpha) |> 
    summarize(pct_greater_d = sum(pct_rel_dist) / n(),
              .groups = "drop")
  
  # Return ----
  return(simulation_tbl_stacked)
}

# Export data ----
c(1,2,4) |>
      # New anonymous function notation
      ## ?`function`
      ### \(x) x + 1
      #### Equivalent to function(x) x + 1 
  map(.f = \(x) ws_simulation(k = x)) |> 
  set_names(nm = c(1, 2, 4)) |> 
  bind_rows(.id = "k") |> 
  mutate(k = as.integer(x = k))  |> 
  write_csv(file = "data/simulation_tbl_stacked.csv")

# Checking ----
read_csv(file = "data/simulation_tbl_stacked.csv")
  