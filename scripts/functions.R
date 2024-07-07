# Libraries ----
library(tidyr)
library(dplyr)
library(purrr)

# Function ws_simulation ----
ws_simulation <- function(mu_1 = 1L,
                          mu_2 = 2L,
                          sigma2_1 = 1L,
                          rho = seq(from = 2L, to = 12L, by = 2L),
                          n_1 = 2L:11L,
                          k = 1L,
                          alpha = c(0.1, 0.05),
                          R = 10000L,
                          seed = 1234L) {
  
  # Parameters ----
  
  ## Means ----
  mu_1 <- mu_1
  mu_2 <- mu_2
  
  ## Variance ----
  sigma2_1 <- sigma2_1
  
  ### Proportions between the variance: sigma2_m_2 / sigma2_1 = rho
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
  null_hypothesis <- mu_1 - mu_2
  
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
  
  ## Seed ----
  set.seed(seed = seed, 
           kind = "Mersenne-Twister", 
           normal.kind = "Inversion", 
           sample.kind = "Rejection")
  
  # Simulations ----
  ## Create samples id: (rho, n_1, r)
  ### Where r = 1, ..., R
  simulation_tbl <- tidyr::expand_grid(rho = rho, 
                                       n_1 = n_1,
                                       r = 1:R)
  
  ## Create samples
  simulation_tbl <- simulation_tbl |> 
    dplyr::mutate(samples_1 = purrr::map2(.x = r, 
                                          .y = n_1,
                                          .f = ~ probdf(n = .y,
                                                        mean = mu_1,
                                                        sd = sqrt(x = sigma2_1))),
                  samples_2 = purrr::pmap(.l = list(r, n_1, rho),
                                          .f = ~ probdf(n = ..2*k,
                                                        mean = mu_2,
                                                        sd = sqrt(x = ..3*sigma2_1))))
  
  ## Calculate mu_1, mu_2, sigma2_1, sigma2_2 for each sample
  simulation_tbl <- simulation_tbl |> 
    dplyr::mutate(samples_mu_1 = purrr::map(.x = samples_1, 
                                            .f = mean),
                  samples_mu_2 = purrr::map(.x = samples_2, 
                                            .f = mean),
                  samples_sigma2_1 = purrr::map(.x = samples_1, 
                                                .f = var),
                  samples_sigma2_2 = purrr::map(.x = samples_2, 
                                                .f = var),
                  .keep = "unused") |> 
    tidyr::unnest(cols = c(samples_mu_1:samples_sigma2_2))
  
  ## Calculate T=t, rho_est
  simulation_tbl <- simulation_tbl |>
    dplyr::mutate(t = t_stat(mu_1 = samples_mu_1,
                             mu_2 = samples_mu_2, 
                             sigma2_1 = samples_sigma2_1, 
                             sigma2_2 = samples_sigma2_2, 
                             n_1 = n_1, n_2 = k*n_1, 
                             null_hypothesis = null_hypothesis)) |> 
    dplyr::mutate(rho_est = samples_sigma2_2/samples_sigma2_1) |> 
    dplyr::select(-c(samples_mu_1:samples_sigma2_2))
  
  ## Calculate rho_est, nu_est, t_alpha/2_n_est (t_ws_cv_lower) 
  simulation_tbl <- rep.int(x = list(simulation_tbl), 
                            times = length(alpha)) |> 
    dplyr::bind_rows() |> 
    dplyr::mutate(alpha = rep(x = alpha, 
                              each = nrow(simulation_tbl))) |>
    dplyr::mutate(nu_est = df_ws(n_1 = n_1,
                                 rho = rho_est,
                                 k = k)) |> 
    dplyr::mutate(t_ws_cv_lower = t_cv_lower(alpha = alpha,
                                             df = nu_est)) |> 
    dplyr::group_by(rho, n_1, alpha) |> 
    dplyr::mutate(alpha_est = ecdf(x = t)(t_ws_cv_lower) + (1 - ecdf(x = t)(-t_ws_cv_lower))) |> 
    dplyr::ungroup()
  
  # Return ----
  return(simulation_tbl)
}

# Function ws_pct_rel_dist ----
ws_pct_rel_dist <- function(simulation_tbl,
                            d = 0.10) {
  
  # Calculate proportion of replicates where 
  # abs(alpha_est - alpha) / alpha > d 
  
  # Parameters ----
  
  ## Lower bound ----
  d <- d
  
  simulation_pct_rel_dist_tbl <- simulation_tbl |> 
    dplyr::mutate(pct_rel_dist = (abs(alpha - alpha_est) / alpha) > d) |> 
    dplyr::group_by( rho, n_1, alpha) |> 
    dplyr::summarize(pct_greater_d = sum(pct_rel_dist) / n(),
                     .groups = "drop")
  
  return(simulation_pct_rel_dist_tbl)
}
