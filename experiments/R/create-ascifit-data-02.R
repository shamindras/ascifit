# Load libraries
# devtools::load_all()
# devtools::document()
# library(fAsianOptions) # needed for the fAsianOptions::erf function
library(tidyverse)

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
gen_asci_data_type1 <- function(n, p, eta, sigma){
  # Rademacher(p) contamination
  xi <- (2 * rbinom(n = n, size = 1, prob = p) - 1)

  # True response
  mu <- seq(from = n * eta, to = n, length.out = n) / n
  eps <- rnorm(n = n, mean = 0, sd = sigma)

  # Contaminated response
  R <- xi * mu + eps

  # ASCI generated data
  asci_gen_data <- tibble::tibble(
    mu_i_eta = mu,
    eps_i = eps,
    xi_i = xi,
    y_i = mu + eps,
    r_i = R,
    t_i = base::abs(R)
  )

  # # Load the ASCI generated data, and eta values into a list
  # out_list <- list("asci_gen_data" = asci_gen_data,
  #                  "eta" = exp_eta)
  base::return(asci_gen_data)
}

ascifit_wrapper <- function(gen_asci_data, eta){
  R_vals <- gen_asci_data$r_i
  out_list <- ascifit(R, eta = eta)
  base::return(out_list)
}

# Test Wrapper ------------------------------------------------------------
# Set seed
set.seed(10747457)
asci_data <- gen_asci_data_type1(n = exp_n,
                                 p = exp_p,
                                 eta = exp_eta,
                                 sigma = exp_sigma)
out_ascifit <- ascifit_wrapper(gen_asci_data = asci_data, eta = exp_eta)
plot(asci_data$r_i)
points(out_ascifit$mu_hat, col = "blue")
points(asci_data$mu_i_eta, col = "red")
dev.off()

# Generate ASCI grid -------------------------------------------------------------
# Define grid parameter value ranges
exp_n <- base::seq(from = 50, to = 100, by = 50)
exp_sigma <- base::seq(from = 0.5, to = 1, by = 0.5)
exp_eta <- 1/5
exp_p <- base::seq(from = 0.25, to = 0.75, by = 0.5)
exp_reps <- 1:5

# Generate the ASCI parameter grid
asci_grid <- tidyr::crossing(n = exp_n,
                             sigma = exp_sigma,
                             eta = exp_eta,
                             p = exp_p,
                             reps = exp_reps) %>%
  dplyr::mutate(sim_no = dplyr::row_number()) %>%
  dplyr::relocate(.data = ., sim_no, reps)

# Unit test: the grid has correct number of rows
length(exp_n) * length(exp_sigma) * length(exp_eta) * length(exp_p) *
  length(exp_reps) - base::nrow(asci_grid)

# Run ASCIFIT over the grid
out1 <- asci_grid %>%
  # dplyr::select(-sim_no, -reps) %>%
  # dplyr::mutate(asci_data = purrr::pmap(.l = ., .f = gen_asci_data_type1),
  dplyr::mutate(asci_data = purrr::pmap(.l = list(.data[["n"]],
                                                  .data[["sigma"]],
                                                  .data[["eta"]],
                                                  .data[["p"]]), .f = gen_asci_data_type1),
                ascifit = purrr::map2(.x = .data[["asci_data"]],
                                      .y = .data[["eta"]],
                                      .f = ~ascifit_wrapper(gen_asci_data = .x,
                                                            eta = .y)),
                mu_hat_inv = purrr::map(.x = .data[["ascifit"]], .f = ~.x$mu_hat_inv),
                sigma_hat = purrr::map(.x = .data[["ascifit"]], .f = ~.x$sigma_hat),
                mu_hat = purrr::map(.x = .data[["ascifit"]], .f = ~.x$mu_hat))
