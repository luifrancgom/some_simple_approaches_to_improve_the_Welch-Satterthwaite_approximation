# Libraries
library(tidyverse)

# Simulations
rho <- seq.int(from = 2, to = 10, by = 2) |> 
  as.integer()
n <- 2:12

expand_grid(rho = rho, n = n)



