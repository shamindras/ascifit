devtools::document()
devtools::load_all()

set.seed(32546)

# Define key functions ---------------------------------------------------------
# TODO: Move these functions into the main R package once we update them with
# roxygen documentation

#' A smooth monotonic univariate \code{sin} function
#' This is based on the paper
#' \url{https://projecteuclid.org/download/pdfview_1/euclid.ejs/1580871776}
#'
#' @param t: (double) a real number on which the function is defined
#'
#' @return (double) function value
#' @export
fnmon_smooth <- function(t){
  base::return(t + base::sin(4 * pi * t)/16)
}

#' Variance of a Folded Normal variable with parameters \eqn{(\mu, \sigma)}
#'
#' Let \eqn{X \sim N(\mu, \sigma)} be an __unfolded__ random variable with mean
#' \eqn{(\mu)} and variance \eqn{(\sigma^{2})}. Then \eqn{T := |X|} is a
#' __Folded Normal__ random variable with parameters \eqn{(\mu, \sigma)}. This
#' function returns the mean of \eqn{T}, i.e.,
#' \deqn{f(\mu, \sigma) := \mathbf{E}(T) = \sigma \sqrt{2 / \pi} \exp \left(-\mu^{2} /\left(2 \sigma^{2}\right)\right)-\mu(1-2 \Phi(\mu / \sigma))}
#'
#' @param mu : The mean parameter of the __unfolded__ Normal variable.
#' @param sigma : The standard of the __unfolded__ Normal variable.
#'
#' @return The mean of the folded normal random variable.
#' @export
f_foldnorm_mean <- function(mu, sigma) {
  out <- sigma * base::sqrt(2 / pi) * base::exp((-1 * mu^2) / (2 * sigma^2)) -
    (mu * (1 - 2 * pnorm(q = (mu / sigma), mean = 0, sd = 1)))
  base::return(out)
}

#' Variance of a Folded Normal variable with parameters \eqn{(\mu, \sigma)}
#'
#' Let \eqn{X \sim N(\mu, \sigma)} be an __unfolded__ random variable with mean
#' \eqn{(\mu)} and variance \eqn{(\sigma^{2})}. Then \eqn{T := |X|} is a
#' __Folded Normal__ random variable with parameters \eqn{(\mu, \sigma)}. This
#' function returns the variance of \eqn{T}, i.e.,
#' \deqn{g(\mu, \sigma) := \mathbf{V}(T) = \mu^{2} + \sigma^{2} - (f(\mu, \sigma))^{2}},
#' where \eqn{f(\mu, \sigma) := \mathbf{E}(T)}, is the mean of the same folded
#' normal variable.
#'
#' @param mu : The mean parameter of the __unfolded__ Normal variable.
#' @param sigma : The standard of the __unfolded__ Normal variable.
#'
#' @return The variance of the folded normal random variable.
#' @export
g_foldnorm_var <- function(mu, sigma) {
  base::return(mu^2 + sigma^2 - (f_foldnorm_mean(mu = mu, sigma = sigma))^2)
}

# Create ascifit simulation data and run ascifit -------------------------------
# Define key ascifit simulation data parameters
n_tot <- 500
norm_mean <- 0
norm_sigma <- 1 # Standard deviation, not variance
eta <- 0.1
p <- 0.5 # sign-flip - Rademacher(p), so +1 with prob p and -1 with prob (1 - p)

# Create ascifit simulation data ----
# This data is normalized on the x-axis over the unit interval for display
# convenience.
sim_df <- base::seq(from = 1, to = n_tot, by = 1) %>%
  tibble::enframe(x = ., name = NULL, value = "n") %>%
  dplyr::mutate(.data = .,
                x_i   = n/n_tot,
                mu_i  = purrr::map_dbl(.x = x_i,
                                       .f = ~fnmon_smooth(t = .x)),
                mu_i_eta = mu_i + eta,
                eps_i = stats::rnorm(n = n_tot,
                                     mean = norm_mean,
                                     sd = norm_sigma),
                xi_i = sample(x = c(-1, 1),
                              size = n_tot,
                              replace = TRUE,
                              prob = c(p, 1 - p)),
                y_i = mu_i_eta + eps_i,
                r_i = xi_i * y_i,
                t_i = base::abs(r_i))

# Run pava on the original data, i.e., on the y_i
# TODO: remove this later when we fit ascifit on the
pava_fit <- stats::isoreg(x = sim_df$x_i, y = sim_df$y_i)

# Estimate (mu, sigma) using ASCIFIT --------------------------------------
# TODO: Write the 3 ascifit steps here

# ASCIFIT Step 1: fit pava on t_i := abs(r_i) terms.
# Get the t_i values
# t_i <- sim_df$t_i

# Fit pava on the t_i values to obtain t_hat_i
# These are the t_hat_i values on the pre-processed data
# t_hat_i_pre <- stats::isoreg(x = sim_df$x_i, y = sim_df$t_i)

# ASCIFIT Step 2: solve sigma satisfying the given expression
# Get the target value of the equation,
# i.e. `\frac{1}{n} \sum_{i = 1}^{n} T_{i}^{2}`
# sig_tgt <- mean((t_i)^2)

# TODO: Add code to solve for sigma_hat
# I believe we need to use stats::uniroot to do this efficiently

# ASCIFIT Step 3: Once sigma_hat is found correct $\hat \absresp_{i}$
# These are the t_hat_i values on the post-processed data
# t_hat_i_post <- TODO: Add code here

# Append the ASCIFIT estimates onto the data
# TODO: The ascifit estimates need to be included here, once we estimate
# them above
sim_df <- sim_df %>%
  dplyr::mutate(.data = .,
                pava_fit_i = pava_fit$yf)

# Plot the simulated results ---------------------------------------------------
sim_plt <- sim_df %>%
  dplyr::select(-eps_i) %>%
  tidyr::pivot_longer(data = .,
                      cols = -c("n", "x_i", "mu_i", "xi_i", "t_i"),
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

# TODO: Checklist
# 1. Check why doesn't isotone::gpava work above?
# pava_fit <- isotone::gpava(z = sim_df$y_i, y = sim_df$x_i, ties = "primary")

# 2. Why does

# Unit tests -------------------------------------------------------------------
# TODO: Move these unit tests into the main R package once we finalize them.

# Tests: f_foldnorm_mean
# mu values
mu_vals <- seq(0, 10, 0.1)

# Unit test 1: # All f(mu, 1) - mu >= 0, fixing sigma = 1
all(f_foldnorm_mean(mu = mu_vals, sigma = 1.5) - mu_vals >= 0)

# Unit test 1: Check f(0, 1) = sqrt(2 / pi)
f_foldnorm_mean(mu = 0, sigma = 1) - sqrt(2 / pi) == 0

# Tests: f_foldnorm_var
g_foldnorm_var(mu = 0, sigma = 1) - (1 - 2 / pi) == 0
