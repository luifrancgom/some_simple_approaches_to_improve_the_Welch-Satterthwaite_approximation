# Libraries ----
library(ggplot2)
library(dplyr)
library(readr)
library(latex2exp)

# Data ----
simulation_tbl <- readr::read_rds(file = "data/simulation_tbl.rds")

## Filter data ----
k <- 1
rho <- 6
n_1 <- 9
alpha <- 0.05

simulation_filter_tbl <- simulation_tbl |> 
  filter(k == !!k, rho == !!rho, n_1 == !!n_1,
         alpha == !!alpha)

# Visualizations ----
## True distribution ----
bins = 55

nu_est_1 <- simulation_filter_tbl$nu_est[1]
t_ws_cv_lower_1 <- simulation_filter_tbl$t_ws_cv_lower[1]
alpha_est_1 <- simulation_filter_tbl$alpha_est[1]

breaks_scale_x <- c(t_ws_cv_lower_1, -t_ws_cv_lower_1)
labels_scale_x <- c(TeX(r'($t_{\widehat{\nu}_1,\frac{\alpha}{2}}=-2.23$)'),
                    TeX(r'($t_{\widehat{\nu}_1,1-\frac{\alpha}{2}}=2.23$)'))

p1 <- simulation_filter_tbl |> 
  ggplot(aes(x = t, 
             y = after_stat(count) / sum(after_stat(count)))) + 
  geom_histogram(bins = bins, 
                 color = "black",
                 closed = "right") + 
  geom_vline(xintercept = breaks_scale_x, 
             color = "#E31A1C") +
  annotate(geom = "text",
           x = 4, y = 0.03,
           label = TeX(r'($\widehat{\alpha}_1=0.0462$)')) +
  annotate(geom = "curve",
           x = 2.7, y = 0.003,
           xend = 4, yend = 0.025,
           arrow = arrow(length = unit(0.1, "inches"),
                         ends = "both"),
           color = "#2C3E50") +
  annotate(geom = "curve",
           x = -2.7, y = 0.003,
           xend = 2.7, yend = 0.03,
           arrow = arrow(length = unit(0.1, "inches"),
                         ends = "both"),
           curvature = -0.3,
           color = "#2C3E50") +
  scale_x_continuous(breaks = breaks_scale_x,
                     labels = labels_scale_x) +
  labs(x = NULL,
       y = "Relative frequency",
       title = TeX(r'($T$ empirical distribution with $R=10000$)'),
       subtitle = TeX(r'($k = 1, \rho = 6, \n_1 = 9, \widehat{\nu}_1=9.945, \alpha=0.05$)'))

p1

## Alpha estimated distribution ----
breaks_scale_x <- c(0.9*alpha, alpha, 1.1*alpha)
labels_scale_x <- c(TeX(r'($0.9\alpha = 0.045$)'),
                    TeX(r'($\alpha = 0.05$)'),
                    TeX(r'($1.1\alpha = 0.055$)'))

p2 <- simulation_filter_tbl |> 
  ggplot(aes(x = alpha_est), 
         y = after_stat(count) / sum(after_stat(count))) + 
  geom_histogram(color = "black",
                 bins = bins) + 
  geom_vline(xintercept = breaks_scale_x[c(1,3)],
             color = "#E31A1C") + 
  geom_vline(xintercept = breaks_scale_x[2],
             color = "#2C3E50")  +
  scale_x_continuous(breaks = breaks_scale_x,
                     labels = labels_scale_x) +
  labs(x = NULL,
       y = "Relative frequency",
       title = TeX(r'($\widehat{\alpha}$ empirical distribution with $R=10000$)'),
       subtitle = TeX(r'($k = 1, \rho = 6, \n_1 = 9, \alpha=0.05$)'))

p2

# Export ----
p1 |> ggsave(filename = "auxiliary_plots/T_empirical_dist.png",
       plot = _,
       width = 7,
       height = 5,
       units = "in")

p2 |> ggsave(filename = "auxiliary_plots/alpha_est_empirical_est.png",
             plot = _,
             width = 7,
             height = 5,
             units = "in")
  
