

library(tidyverse)

sd <- 1 # sigma
mu <- 10

# mean_foldn <- function(mu, sd) {
#   out <- mu^2 + sd^2 - ((sqrt(2 / pi)) * (sd) * exp((-1 * mu^2) / (2*sd^2)) +
#     mu * (1 - 2 * pnorm((-mu) / sd)))^2
#   return(out)
# }

mean_foldn <- function(mu, sd) {
  out <- sqrt(2 / pi) * (sd) * exp((-1 * mu^2) / (2*sd^2)) -
    mu * (1 - 2 * pnorm((mu) / sd))
  return(out)
}

test_out <- tidyr::crossing(
  mu = seq(from = 0, to = 5, by = 0.1),
  sd = seq(from = 0.1, to = 10, by = 0.5)
) %>%
  dplyr::mutate(
    mean_foldn = purrr::pmap_dbl(
      .l = .,
      .f = mean_foldn
    ),
    sd_fac = forcats::as_factor(x = as.character(sd))
  ) %>%
  dplyr::relocate(.data = ., sd_fac, .before = mean_foldn)
test_out

plt_out <- test_out %>%
  ggplot2::ggplot(data = .,
                  mapping = aes(x = mu,
                                y = mean_foldn,
                                col = sd_fac)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  NULL
plt_out
plotly::ggplotly(plt_out)


library(tidyverse)

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
