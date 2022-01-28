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
  # R <- xi * mu + eps
  R <- xi * (mu + eps)

  # Fit ASCIFIT to the data
  out_ascifit <- ascifit(R_vals = R, eta = eta)

  # ASCI generated data
  asci_gen_data <- tibble::tibble(
    x_i = seq(from = n * eta, to = n, length.out = n),
    mu_i = mu,
    eps_i = eps,
    xi_i = xi,
    y_i = mu + eps,
    r_i = R,
    t_i = base::abs(R),
    mu_hat_inv_i = out_ascifit$mu_hat_inv,
    mu_hat_i = out_ascifit$mu_hat
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

mse <- function(estimate, actual){
  base::return(mean((estimate - actual)**2))
}

# Generate ASCI grid -------------------------------------------------------------
# Define grid parameter value ranges
exp_n <- base::seq(from = 50, to = 200, by = 50)
exp_sigma <- base::seq(from = 0.5, to = 1, by = 0.5)
exp_eta <- 1/5
exp_p <- base::seq(from = 0.25, to = 0.75, by = 0.5)
exp_reps <- 1:5
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
tot_sims - base::nrow(asci_grid)

# Run ASCIFIT over the grid
out_ascifit1 <- asci_grid %>%
  dplyr::mutate(asci_data = purrr::pmap(.l = list(.data[["n"]],
                                                  .data[["sigma"]],
                                                  .data[["eta"]],
                                                  .data[["p"]]),
                                        .f = gen_asci_data_type1))
out_ascifit2 <- out_ascifit1 %>%
  dplyr::mutate(mu_hat_inv = purrr::map(.x = .data[["asci_data"]], .f = ~.x$mu_hat_inv_i),
                sigma_hat = purrr::map(.x = .data[["asci_data"]], .f = ~.x$sigma_hat_i),
                mu_hat = purrr::map(.x = .data[["asci_data"]], .f = ~.x$mu_hat_i),
                mu = purrr::map(.x = .data[["asci_data"]], .f = ~.x$mu_i),
                mse = purrr::map2_dbl(.x = .data[["mu_hat"]], .y = .data[["mu"]],
                                      .f = ~mse(estimate = .x, actual = .y)))
out_ascifit2$mse


# Plot mean-MSE of ASCIFIT over grid --------------------------------------
# Wrangle the plotting data
# For a fixed value of eta, p, summarize the mean-MSE by n and sigma
out_ascifit2_mse <- out_ascifit2 %>%
  dplyr::select(n, eta, p, sigma, mse) %>%
  # dplyr::filter(eta == 1/5, p == 0.25, sigma == 1) %>%
  dplyr::filter(eta == 1/5, p == 0.25) %>%
  group_by(n, sigma) %>%
  dplyr::summarize(mean_mse = base::mean(mse),
                   sd_mse = stats::sd(mse)) %>%
  dplyr::mutate(se_mse = sd_mse/base::sqrt(n),
                sigma = as.factor(sigma),
                n = as.factor(n))

# For a fixed value of eta, p, plot the mean-MSE by n and sigma

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- ggplot2::position_dodge(0.1) # move them .05 to the left and right

out_ascifit2_mse %>%
  ggplot2::ggplot(data = ., ggplot2::aes(x = n, y = mean_mse, color = sigma)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_mse - 2*se_mse,
                                      ymax = mean_mse + 2*se_mse), width = .1,
                         position = pd) +
  ggplot2::geom_line(position=pd) +
  ggplot2::geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  ggplot2::labs(title = "Mean MSE for estimating mu using ASCIFIT",
                x = "Sample size (n)",
                y = latex2exp::TeX("Sample Mean-MSE $\\| \\hat{\\mu} - mu \\|^{2}$")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.justification=c(1,0),
                 legend.position=c(1,0))               # Position legend in bottom right

out_ascifit3 <- out_ascifit2 %>%
  dplyr::filter(eta == 1/5, p == 0.25, sigma == 1, reps == 1, n == 200)

out_ascifit3$asci_data
out_ascifit3$ascifit
out_ascifit2

# Plot the simulated results ---------------------------------------------------
sim_plt <- out_ascifit3 %>%
  dplyr::pull(asci_data) %>%
  .[[1]] %>%
  dplyr::select(-eps_i) %>%
  tidyr::pivot_longer(data = .,
                      cols = -c("mu_i", "xi_i", "t_i"),
                      names_to = c("y_type")) %>%
  dplyr::mutate(y_type = base::as.factor(y_type)) %>%
  ggplot2::ggplot(data = .,
                  mapping = ggplot2::aes(x = x_i, y = value)) +
  # ggplot2::geom_line(ggplot2::aes(color = y_type)) +
  ggplot2::geom_point(ggplot2::aes(color = y_type)) +
  ggplot2::scale_color_manual(values = c("purple", "blue", "orange", "darkgreen"))  +
  ggplot2::ylim(-5, 5) +
  ggplot2::xlim(0, 1) +
  ggplot2::labs(
    title = "Plot of true values, response, and fitted PAVA",
    x = "x",
    y = "true values, response, and fitted PAVA"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.title = ggplot2::element_blank())
sim_plt

