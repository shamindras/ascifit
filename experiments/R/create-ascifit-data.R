devtools::document()
devtools::load_all()

set.seed(32546)

# Define key functions ---------------------------------------------------------
#' A smooth monotonic univariate \code{sin} function
#' This is based on the paper \url{https://projecteuclid.org/download/pdfview_1/euclid.ejs/1580871776}
#'
#' @param t: (double) a real number on which the function is defined
#'
#' @return (double) function value
#' @export
fnmon_smooth <- function(t){
  base::return(t + base::sin(4 * pi * t)/16)
}

n_tot <- 500
norm_mean <- 0
norm_sigma <- 1 # Standard deviation, not variance
eta <- 0.1
p <- 0.5 # sign-flip - Rademacher(p), so +1 with prob p and -1 with prob (1 - p)

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
                r_i = xi_i * y_i)

pava_fit <- stats::isoreg(x = sim_df$x_i, y = sim_df$y_i)

# TODO: Write the 3 ascifit steps here
# asci_fit_step1 <- stats::isoreg(x = sim_df$x_i, y = base::abs(sim_df$r_i))

sim_df   <- sim_df %>%
  dplyr::mutate(.data = .,
                pava_fit_i = pava_fit$yf)

sim_plt <- sim_df %>%
  dplyr::select(-eps_i) %>%
  tidyr::pivot_longer(data = .,
                      cols = -c("n", "x_i", "mu_i", "xi_i"),
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


# TODO: Check why doesn't isotone::gpava work above?
# pava_fit <- isotone::gpava(z = sim_df$y_i, y = sim_df$x_i, ties = "primary")
