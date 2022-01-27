# Set seed
set.seed(10747457)

# Setup grid parameters
exp_n <- 50
exp_p <- 0.5
exp_eta <- 1/5
exp_sigma <- 1

# ASCI generated data -----------------------------------------------------
# Rademacher(p) contamination
xi <- (2 * rbinom(n = exp_n, size = 1, prob = exp_p) - 1)

# True response
mu <- seq(from = exp_n * exp_eta, to = exp_n, length.out = exp_n) / exp_n
eps <- rnorm(n = exp_n, mean = 0, sd = exp_sigma)

# Contaminated response
R <- xi * mu + eps

# Check that the vectors are of the same length
length(mu)
length(R)

# run ASCIFIT on generated data -------------------------------------------
test <- ascifit(R, eta = exp_eta)
length(test$mu_hat)

# Plot results ------------------------------------------------------------
# dev.off()
plot(R)
points(test$mu_hat, col = "blue")
points(mu, col = "red")
dev.off()

# Calculate MSE -----------------------------------------------------------
mean((mu - test$mu_hat)**2)
exp_sigma - test$sigma_hat
