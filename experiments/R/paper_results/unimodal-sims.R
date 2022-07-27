# Load libraries
# devtools::load_all()
# devtools::document()
# library(fAsianOptions) # needed for the fAsianOptions::erf function
library(tidyverse)
library(here)
library(scales)
library(ggthemes)

# Function to produce MSE
mse <- function(estimate, actual) {
  base::return(mean((estimate - actual)**2))
}

# Need to write the erf explicitly, since we can't load fAsianOptions::erf
# on R version 4.2.1+
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

inverse_mean_sigma <- function(f, sigma, lower = 0, upper = 10000) {
  function(y) uniroot((function(x) f(x, sigma) - y), lower = lower, upper = upper)[1]
}

# this function is the same as the one below
mean_theta_sigma <- function(theta, sigma) {
  x <- ((theta / sigma)**2 / 2)
  base::return((sqrt(pi) * sqrt(x) * erf(sqrt(x)) + exp(-x)) * sqrt(2 / pi) * sigma)
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
  out_list <- list(
    "mu_hat_inv" = theta1,
    "mu_hat" = theta2,
    "sigma_hat" = sigma
  )
  base::return(out_list)
}

# Wrapper - generate ASCI data --------------------------------------------
gen_asci_data_type1 <- function(tot_sims, sim_no, n, p, eta, sigma) {
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

ascifit_wrapper <- function(gen_asci_data, eta) {
  R_vals <- gen_asci_data$r_i
  out_list <- ascifit(R_vals = R_vals, eta = eta)
  base::return(out_list)
}

asci_data_wrapper <- function(asci_data_out, ascifit_out) {
  out_comb_data <- asci_data_out %>%
    dplyr::mutate(
      mu_hat_inv_i = ascifit_out$mu_hat_inv,
      mu_hat_i = ascifit_out$mu_hat
    )
  base::return(out_comb_data)
}

# Set seed for reproducibility
set.seed(551243)

# Generate ASCI grid -------------------------------------------------------------
# Define grid parameter value ranges
exp_n <- c(1000L)
exp_sigma <- c(1.5)
exp_eta <- 1 / 5
exp_p <- 0.5
exp_reps <- 1
# Get total number of sims. Useful for printing
tot_sims <- length(exp_n) * length(exp_sigma) * length(exp_eta) * length(exp_p) * length(exp_reps)

# Generate the ASCI parameter grid
asci_grid <- tidyr::crossing(
  n = exp_n,
  sigma = exp_sigma,
  eta = exp_eta,
  p = exp_p,
  reps = exp_reps
) %>%
  # Add in row number index for the replications
  dplyr::mutate(
    sim_no = dplyr::row_number(),
    tot_sims = tot_sims
  ) %>%
  dplyr::relocate(.data = ., tot_sims, sim_no, reps)

# Unit test: the grid has correct number of rows
assertthat::assert_that(tot_sims - base::nrow(asci_grid) == 0)

# Run ASCIFIT over the grid
out_ascifit1 <- asci_grid %>%
  dplyr::mutate(
    asci_data_out = purrr::pmap(
      .l = list(
        .data[["tot_sims"]],
        .data[["sim_no"]],
        .data[["n"]],
        .data[["p"]],
        .data[["eta"]],
        .data[["sigma"]]
      ),
      .f = gen_asci_data_type1
    ),
    ascifit_out = purrr::map2(
      .x = .data[["asci_data_out"]],
      .y = .data[["eta"]],
      .f = ~ ascifit_wrapper(
        gen_asci_data = .x,
        eta = .y
      )
    )
  ) %>%
  dplyr::mutate(ascifit_comb = purrr::pmap(
    .l = list(
      .data[["asci_data_out"]],
      .data[["ascifit_out"]]
    ),
    .f = asci_data_wrapper
  ))

# Produce plots
# Produce plots for conference paper
out_ascifit2 <- out_ascifit1 %>%
  dplyr::mutate(
    mu_hat_inv = purrr::map(.x = .data[["ascifit_out"]], .f = ~ .x$mu_hat_inv),
    sigma_hat = purrr::map(.x = .data[["ascifit_out"]], .f = ~ .x$sigma_hat),
    mu_hat = purrr::map(.x = .data[["ascifit_out"]], .f = ~ .x$mu_hat),
    mu = purrr::map(.x = .data[["asci_data_out"]], .f = ~ .x$mu_i),
    mse = purrr::map2_dbl(
      .x = .data[["mu_hat"]], .y = .data[["mu"]],
      .f = ~ mse(estimate = .x, actual = .y)
    )
  )


# Plot the simulated results ---------------------------------------------------
plt_eta <- 1 / 5
plt_p <- 0.5
plt_sigma <- 1.5
plt_reps <- 1
plt_n <- 1000
out_ascifit3 <- out_ascifit2 %>%
  dplyr::filter(
    eta == plt_eta, p == plt_p, sigma == plt_sigma,
    reps == plt_reps, n == plt_n
  )
out_ascifit3

out_ascifit3 %>%
  dplyr::pull(ascifit_comb) %>%
  .[[1]]

sim_plt <- out_ascifit3 %>%
  dplyr::pull(ascifit_comb) %>%
  .[[1]] %>%
  dplyr::select(-eps_i) %>% # -mu_hat_inv_i
  tidyr::pivot_longer(
    data = .,
    cols = -c("x_i", "xi_i", "t_i"),
    names_to = c("y_type")
  ) %>%
  dplyr::mutate(y_type = base::as.factor(y_type)) %>%
  ggplot2::ggplot(
    data = .,
    mapping = ggplot2::aes(x = x_i, y = value)
  ) +
  ggplot2::geom_point(ggplot2::aes(color = y_type)) +
  ggplot2::ylim(-5, 5) +
  ggplot2::xlim(0, plt_n) +
  # ggthemes::theme_economist() +
  ggthemes::theme_igray() +
  ggplot2::scale_color_hue(labels = c(
    latex2exp::TeX("$\\hat{\\mu}$"),
    latex2exp::TeX("$\\hat{\\mu}_{naive}$"),
    latex2exp::TeX("$\\mu$"),
    latex2exp::TeX("$R$"),
    latex2exp::TeX("$Y$")
  )) +
  ggplot2::labs(
    title = "Plot of true values, response, and ASCIFIT",
    x = "x",
    y = "response"
    # y = "true values, response, and ASCIFIT"
  ) +
  # ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom",
    text = element_text(size = 20),
    legend.title = ggplot2::element_blank()
  )
sim_plt
