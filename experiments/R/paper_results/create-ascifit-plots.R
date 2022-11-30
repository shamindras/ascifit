# Clean up console and dev list ------------------------------------------------
# clean console
cat("\014")

# clean dev.list, i.e., any pre-generated plots
if(!is.null(dev.list())){
  dev.off()
}

# Load libraries ---------------------------------------------------------------
library(tidyverse)
# library(here)
# library(scales)
# library(ggthemes)

# Define key functions ---------------------------------------------------------
# Compute the MSE between actual vs. estimated values
mse <- function(estimate, actual){
  base::return(mean((estimate - actual)**2))
}

# Produce a more lightweight TikZ representation of plot
gen_tikz_plot <- function(width, height, plt, plt_outdir, plt_outname){
  # turn of any plot devices, so that we save a clean figure output
  if(!is.null(dev.list())){
    dev.off()
  }

  plt_outpath <- here::here(plt_outdir, glue::glue("{plt_outname}.tex"))
  # Delete the TeX file, so that we always recreate it
  if(file.exists(plt_outpath)){
    file.remove(plt_outpath)
  }
  tikzDevice::tikz(
    file = plt_outpath,
    width = width,
    height = width
  )
  print(plt)
  dev.off()

  # removes unnecessary whitespace around the plot when importing into LaTeX
  # this deletes all lines that invisibly mess up the bounding box
  # source: https://stackoverflow.com/a/41186942/4687531
  lines <- readLines(con = plt_outpath)
  lines <- lines[-which(grepl("\\path\\[clip\\]*", lines,
                              perl = F))]
  lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines,
                              perl = F))]
  writeLines(lines, con = plt_outpath)
}

# Set seed for reproducibility -------------------------------------------------
set.seed(523525)

# Define global variables ------------------------------------------------------
# directory where simulated data is stored
DAT_DIR <- here::here("experiments", "R", "paper_results")

# directory to save plots
PLOT_DIR <- here::here("experiments", "R", "paper_results", "plots")

# Data processing for (simulation) figures -------------------------------------
# Read in generated (simulation) ascifit data
out_ascifit1_r <- readr::read_rds(
  file = here::here(DAT_DIR, "out_ascifit1_50_reps.rds")
)

# Clean the simulated data for plotting
out_ascifit2 <- out_ascifit1_r %>%
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

# Plot mean-MSE of ASCIFIT over grid -------------------------------------------
# Wrangle the plotting data
# For a fixed value of eta, p, summarize the mean-MSE by n and sigma
out_ascifit2_mse <- out_ascifit2 %>%
  dplyr::select(n, eta, p, sigma, mse) %>%
  # dplyr::filter(eta == 1/5, p == 0.25, sigma == 1) %>%
  dplyr::filter(eta == 1 / 5, p == 0.5) %>%
  group_by(n, sigma) %>%
  dplyr::summarize(
    mean_mse = base::mean(mse),
    sd_mse = stats::sd(mse)
  ) %>%
  dplyr::mutate(
    se_mse = sd_mse / base::sqrt(n),
    sigma = as.factor(sigma),
    n = as.factor(n)
  )

# For a fixed value of eta, p, plot the mean-MSE by n and sigma

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- ggplot2::position_dodge(0.1) # move them .05 to the left and right

out_ascifit2_mse_plt <- out_ascifit2_mse %>%
  ggplot2::ggplot(data = ., ggplot2::aes(
    x = n, y = mean_mse, color = sigma,
    fill = sigma
  )) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_mse - 2 * se_mse,
      ymax = mean_mse + 2 * se_mse
    ),
    width = .5,
    position = pd,
    linewidth = 0.8
  ) +
  ggplot2::geom_line(position = pd) +
  ggplot2::geom_point(position = pd, size = 2, shape = 21) + # 21 is filled circle
  ggplot2::labs(
    # title = latex2exp::TeX("Sample mean-MSE for estimating $\\mu$ using ASCIFIT"),
    # x = latex2exp::TeX("Sample size $(n)$"),
    # y = latex2exp::TeX("mean-MSE $\\frac{1}{n} \\| \\hat{\\mu} - mu \\|^{2}$")
    # title = "Sample mean-MSE for estimating $\\mu$ using ASCIFIT",
    x = "Sample size $(n)$",
    y = "$\\frac{1}{n} \\| \\hat{\\mu} - \\mu \\|^{2}$",
  ) +
  ggplot2::scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  # ggthemes::theme_economist() +
  # ggthemes::theme_igray() +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  # ggplot2::guides(fill=guide_legend(title="\\sigma")) +
  NULL

# Inspect plot
out_ascifit2_mse_plt

# Produce a png of the plot
cowplot::ggsave2(filename = here::here(PLOT_DIR, "out_ascifit1_50_reps.png"),
                 plot = out_ascifit2_mse_plt,
                 # width = 10,
                 # height = 50,
                 dpi = 900)

# Produce a more lightweight TikZ representation of plot
gen_tikz_plot(width = 5,
              height = 8,
              plt = out_ascifit2_mse_plt,
              plt_outdir = PLOT_DIR,
              plt_outname = "out_ascifit1_50_reps")

# Plot the simulated results ---------------------------------------------------
plt_eta <- 0.2
plt_p <- 0.5
plt_sigma <- 1.5
plt_reps <- 33
plt_n <- 1000
out_ascifit3 <- out_ascifit2 %>%
  dplyr::filter(eta == plt_eta, p == plt_p, sigma == plt_sigma,
                reps == plt_reps, n == plt_n)
out_ascifit3

out_ascifit3 %>%
  dplyr::pull(ascifit_comb) %>%
  .[[1]]

sim_plt <- out_ascifit3 %>%
  dplyr::pull(ascifit_comb) %>%
  .[[1]] %>%
  dplyr::select(-eps_i) %>% # -mu_hat_inv_i
  tidyr::pivot_longer(data = .,
                      cols = -c("x_i", "xi_i", "t_i"),
                      names_to = c("y_type")) %>%
  dplyr::mutate(y_type = base::as.factor(y_type)) %>%
  ggplot2::ggplot(data = .,
                  mapping = ggplot2::aes(x = x_i, y = value)) +
  ggplot2::geom_point(ggplot2::aes(color = y_type)) +
  ggplot2::ylim(-5, 5) +
  ggplot2::xlim(0, plt_n) +
  # ggthemes::theme_economist() +
  # ggthemes::theme_igray() +
  ggplot2::theme_bw() +
  # ggplot2::scale_color_hue(labels = c(latex2exp::TeX("$\\hat{\\mu}$"),
  #                                     latex2exp::TeX("$\\hat{\\mu}_{naive}$"),
  #                                     latex2exp::TeX("$\\mu$"),
  #                                     latex2exp::TeX("$R$"),
  #                                     latex2exp::TeX("$Y$"))) +
  ggplot2::scale_color_hue(labels = c("$\\hat{\\mu}_{\\textnormal{ascifit}}$",
                                      "$\\hat{\\mu}_{\\textnormal{naive}}$",
                                      "$\\mu$",
                                      "$R$",
                                      "$Y$")) +
  ggplot2::labs(
    # title = "Plot of true values, response, and ASCIFIT",
    # x = "$x$",
    x = "",
    y = "response"
    # y = "true values, response, and ASCIFIT"
  ) +
  # ggplot2::theme_minimal() +
  ggplot2::theme(legend.position="bottom",
                 text = element_text(size=12),
                 legend.title = ggplot2::element_blank())
sim_plt

# Produce a png of the plot
cowplot::ggsave2(filename = here::here(PLOT_DIR, "sim_plt_50_reps.png"),
                 plot = sim_plt,
                 # width = 10,
                 # height = 50,
                 dpi = 900)

# Produce a more lightweight TikZ representation of plot
gen_tikz_plot(width = 5,
              height = 8,
              plt = sim_plt,
              plt_outdir = PLOT_DIR,
              plt_outname = "sim_plt_50_reps")
