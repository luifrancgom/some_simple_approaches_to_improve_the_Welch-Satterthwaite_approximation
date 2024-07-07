# Libraries ----
library(purrr)
library(dplyr)
library(readr)

# Functions ----
source(file = "scripts/functions.R")

# Export data ----
simulation_lst <- c(1L, 2L, 4L) |>
  purrr::map(.f = \(x) ws_simulation(k = x)) 

simulation_tbl <- simulation_lst |> 
  dplyr::bind_rows(.id = "k") |> 
  dplyr::mutate(k = as.integer(x = k))

simulation_pct_rel_dist_tbl <- simulation_lst |> 
  purrr::map(.f = \(x) ws_pct_rel_dist(simulation_tbl = x))  |> 
  setNames(nm = c(1L, 2L, 4L)) |> 
  dplyr::bind_rows(.id = "k") |> 
  dplyr::mutate(k = as.integer(x = k))

## Zip file ----

simulation_pct_rel_dist_tbl |> 
  readr::write_csv(file = "data/simulation_pct_rel_dist_tbl.csv")

# Create a connection to a gzipped file
simulation_tbl |>  
  readr::write_rds(file = "data/simulation_tbl.rds",
                   compress = "gz")

# Checking ----
readr::read_csv(file = "data/simulation_pct_rel_dist_tbl.csv")

readr::read_rds(file = "data/simulation_tbl.rds")
