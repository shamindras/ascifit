# Load libraries
# devtools::load_all()
# devtools::document()
library(fAsianOptions) # needed for the fAsianOptions::erf function
library(tidyverse)
library(here)
library(scales)

inverse_mean_sigma <- function (f, sigma, lower = 0, upper = 10000) {
  function (y) uniroot((function (x) f(x, sigma) - y), lower = lower, upper = upper)[1]
}

# this function is the same as the one below
mean_theta_sigma <- function(theta, sigma) {
  x <- ((theta / sigma)**2 / 2)
  base::return((sqrt(pi) * sqrt(x) * fAsianOptions::erf(sqrt(x)) + exp(-x)) * sqrt(2 / pi) * sigma)
}

ascifit <- function(R_vals, eta) {
  # ASCIFIT Step 1: Pre-processing and PAVA
  theta1 <- isoreg(abs(R_vals))$yf
  M <- mean(R_vals**2)

  # ASCIFIT Step 2: Second moment matching
  # Binary search solve for sigma_hat
  L <- 0.001
  R <- sqrt(M)
  subtract <- Inf
  sigma <- (L + R) / 2

  while (abs(M - subtract) > 1e-5) {
    theta1.maxed.eta <- sapply(theta1, function(el) max(el, mean_theta_sigma(eta, sigma)))
    inv <- Vectorize(inverse_mean_sigma(mean_theta_sigma, sigma = sigma, lower = 0))
    subtract <- mean(unlist(inv(theta1.maxed.eta))**2) + sigma**2

    if (subtract > M) {
      R <- sigma
    } else if (subtract < M) {
      L <- sigma
    }
    sigma <- (L + R) / 2
  }

  # ASCIFIT Step 3: Post-processing via plug-in
  # Use sigma_hat plug_in and invert for mu_hat
  theta1.maxed.eta <- sapply(theta1, function(el) max(el, mean_theta_sigma(eta, sigma)))
  inv <- Vectorize(inverse_mean_sigma(mean_theta_sigma, sigma = sigma, lower = 0))
  theta2 <- unlist(inv(theta1.maxed.eta))
  out_list <- list("mu_hat_inv" = theta1,
                   "mu_hat" = theta2,
                   "sigma_hat" = sigma)
  base::return(out_list)
}

# Wrapper - generate ASCI data --------------------------------------------
gen_asci_data_type1 <- function(tot_sims, sim_no, n, p, eta, sigma){
  print(glue::glue("Simulation {sim_no} of {tot_sims}: Progress {scales::percent(sim_no/tot_sims)}\\n"))
  # Rademacher(p) contamination
  xi <- 2 * rbinom(n = n, size = 1, prob = p) - 1

  # True response
  mu <- seq(from = n * eta, to = n, length.out = n) / n
  eps <- rnorm(n = n, mean = 0, sd = sigma)

  # Contaminated response
  # R <- xi * mu + eps
  y <- mu + eps
  R <- xi * y
  T_vals <- base::abs(R)

  # ASCI generated data
  asci_gen_data <- tibble::tibble(
    x_i = seq(from = 1, to = n, length.out = n),
    mu_i = mu,
    eps_i = eps,
    xi_i = xi,
    y_i = y,
    r_i = R,
    t_i = T_vals
  )

  # # Load the ASCI generated data, and eta values into a list
  # out_list <- list("asci_gen_data" = asci_gen_data,
  #                  "eta" = exp_eta)
  base::return(asci_gen_data)
}

ascifit_wrapper <- function(gen_asci_data, eta){
  R_vals <- gen_asci_data$r_i
  out_list <- ascifit(R_vals = R_vals, eta = eta)
  base::return(out_list)
}

asci_data_wrapper <- function(asci_data_out, ascifit_out){
  out_comb_data <- asci_data_out %>%
    dplyr::mutate(mu_hat_inv_i = ascifit_out$mu_hat_inv,
                  mu_hat_i = ascifit_out$mu_hat)
  base::return(out_comb_data)
}

# Set seed for reproducibility
set.seed(4265436)

# Generate ASCI grid -------------------------------------------------------------
# Define grid parameter value ranges
exp_n <- c(100L, 250L, 500L, 1000L)
exp_sigma <- base::seq(from = 0.5, to = 2, by = 0.5)
exp_eta <- 1/5
exp_p <- 0.5
exp_reps <- 1:50
# Get total number of sims. Useful for printing
tot_sims <- length(exp_n) * length(exp_sigma) * length(exp_eta) * length(exp_p) * length(exp_reps)

# Generate the ASCI parameter grid
asci_grid <- tidyr::crossing(n = exp_n,
                             sigma = exp_sigma,
                             eta = exp_eta,
                             p = exp_p,
                             reps = exp_reps) %>%
  # Add in row number index for the replications
  dplyr::mutate(sim_no = dplyr::row_number(),
                tot_sims = tot_sims) %>%
  dplyr::relocate(.data = ., tot_sims, sim_no, reps)

# Unit test: the grid has correct number of rows
assertthat::assert_that(tot_sims - base::nrow(asci_grid) == 0)

# Run ASCIFIT over the grid
out_ascifit1 <- asci_grid %>%
  dplyr::mutate(asci_data_out = purrr::pmap(.l = list(.data[["tot_sims"]],
                                                      .data[["sim_no"]],
                                                      .data[["n"]],
                                                      .data[["p"]],
                                                      .data[["eta"]],
                                                      .data[["sigma"]]),
                                            .f = gen_asci_data_type1),
                ascifit_out = purrr::map2(.x = .data[["asci_data_out"]],
                                      .y = .data[["eta"]],
                                      .f = ~ascifit_wrapper(gen_asci_data = .x,
                                                            eta = .y))) %>%
  dplyr::mutate(ascifit_comb = purrr::pmap(.l = list(.data[["asci_data_out"]],
                                                     .data[["ascifit_out"]]),
                                           .f = asci_data_wrapper))

# Write out the tibble to local directory
readr::write_rds(x = out_ascifit1,
                 file = here::here("experiments", "R", "out_ascifit1_50_reps.rds"),
                 compress = "xz", compression = 9L)
