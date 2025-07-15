# implemented Guassian Process model
# Adapted from here:
# https://peterroelants.github.io/posts/gaussian-process-tutorial/


# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP1
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# Define parameters
total_time <- 1
nb_steps <- 75
delta_t <- total_time / nb_steps
nb_processes <- 5
mean <- 0
stdev <- sqrt(delta_t)

# Simulate the Brownian motions
# Each row will be one process (like shape (nb_processes, nb_steps) in numpy)
steps_matrix <- matrix(
  rnorm(nb_processes * nb_steps, mean, stdev),
  nrow = nb_processes, ncol = nb_steps
)

# Compute cumulative sums along time (columns)
distances <- t(apply(steps_matrix, 1, cumsum))

# Time axis
t <- seq(0, total_time - delta_t, by = delta_t)

# Plot
plot(
  t, distances[1, ], type = "l", col = 1,
  xlim = c(0, 1), ylim = range(distances),
  xlab = "t (time)", ylab = "d (position)",
  main = "Brownian motion process\nPosition over time for 5 independent realizations"
)
for (i in 2:nb_processes) {
  lines(t, distances[i, ], col = i)
}


# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP 2
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

exponentiated_quadratic <- function(xa, xb) {
  # Squared Euclidean distances
  dists <- as.matrix(dist(rbind(xa, xb)))^2
  n_a <- nrow(as.matrix(xa))
  n_b <- nrow(as.matrix(xb))
  sq_dists <- dists[1:n_a, (n_a + 1):(n_a + n_b)]
  return( exp(-0.5 * sq_dists) )
}

# Set limits and grid
xlim <- c(-3, 3)
X <- matrix(seq(xlim[1], xlim[2], length.out = 25), ncol = 1)

# Plot covariance matrix
image(
  1:25, 1:25, Sigma,
  col = heat.colors(100),
  main = "Exponentiated quadratic\nexample of covariance matrix",
  xlab = "x", ylab = "x",
  axes = FALSE
)
ticks <- seq(xlim[1], xlim[2], by = 1)
axis(1, at = seq(1, 25, length.out = length(ticks)), labels = ticks)
axis(2, at = seq(1, 25, length.out = length(ticks)), labels = ticks)
box()

# Covariance matrix
Sigma <- exponentiated_quadratic(X, X)

xlim2 <- c(-4, 4)
X2 <- matrix(seq(xlim2[1], xlim2[2], length.out = 50), ncol = 1)
zero <- matrix(0, ncol = 1)
Sigma0 <- exponentiated_quadratic(X2, zero)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Heatmap
image(
  1:25, 1:25, Sigma,
  col = heat.colors(100),
  main = "Exponentiated quadratic\nexample of covariance matrix",
  xlab = "x", ylab = "x",
  axes = FALSE
)
axis(1, at = seq(1, 25, length.out = length(ticks)), labels = ticks)
axis(2, at = seq(1, 25, length.out = length(ticks)), labels = ticks)
box()

# Covariance curve
plot(
  X2, Sigma0,
  type = "l", lwd = 2,
  xlab = "x", ylab = "covariance",
  main = "Exponentiated quadratic covariance\nbetween x and 0",
  xlim = xlim2
)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP 2
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////


# Exponentiated quadratic kernel
exponentiated_quadratic <- function(xa, xb) {
  dists <- as.matrix(dist(rbind(xa, xb)))^2
  n_a <- nrow(as.matrix(xa))
  n_b <- nrow(as.matrix(xb))
  sq_dists <- dists[1:n_a, (n_a + 1):(n_a + n_b)]
  exp(-0.5 * sq_dists)
}

# Parameters
nb_of_samples <- 41
number_of_functions <- 5
X <- matrix(seq(-4, 4, length.out = nb_of_samples), ncol = 1)
Sigma <- exponentiated_quadratic(X, X)

max(abs(Sigma - t(Sigma)))
eigen(Sigma)$values

Sigma <- Sigma + diag(1e-8, nrow(Sigma))

# Sample from MVN(0, Sigma)
set.seed(123)
L <- chol(Sigma)
ys <- matrix(NA, nrow = number_of_functions, ncol = nb_of_samples)
for (i in 1:number_of_functions) {
  z <- rnorm(nb_of_samples)
  ys[i, ] <- t(L) %*% z
}

# Plot
plot(
  X, ys[1, ],
  type = "o", pch = 16, cex = 0.5, col = 1,
  xlim = c(-4, 4),
  ylim = range(ys),
  xlab = "x",
  ylab = "y = f(x)",
  main = "5 different function realizations at 41 points\nsampled from a GP with exponentiated quadratic kernel"
)
for (i in 2:number_of_functions) {
  lines(X, ys[i, ], type = "o", pch = 16, cex = 0.5, col = i)
}

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP 3
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////


generate_surface <- function(mean, covariance, surface_resolution) {
  x1s <- seq(-5, 5, length.out = surface_resolution)
  x2s <- seq(-5, 5, length.out = surface_resolution)
  pdf <- matrix(0, nrow = surface_resolution, ncol = surface_resolution)

  # In base R, dmvnorm is not built-in.
  # We'll implement directly using density formula
  inv_cov <- solve(covariance)
  det_cov <- det(covariance)
  norm_const <- 1 / (2 * pi * sqrt(det_cov))

  for (i in 1:surface_resolution) {
    for (j in 1:surface_resolution) {
      x <- c(x1s[i], x2s[j])
      pdf[i, j] <- norm_const * exp(-0.5 * t(x - mean) %*% inv_cov %*% (x - mean))
    }
  }

  list(x1 = x1s, x2 = x2s, pdf = pdf)
}

X_strong <- matrix(c(0, 0.2), ncol = 1)
Σ_strong <- exponentiated_quadratic(X_strong, X_strong)

X_weak <- matrix(c(0, 2), ncol = 1)
Σ_weak <- exponentiated_quadratic(X_weak, X_weak)

mean_ <- c(0, 0)
surface_resolution <- 50

surf_strong <- generate_surface(mean_, Σ_strong, surface_resolution)
surf_weak   <- generate_surface(mean_, Σ_weak, surface_resolution)

X_0_index  <- which(abs(X - 0) < 1e-8)
X_02_index <- which(abs(X - 0.2) < 1e-8)
X_2_index  <- which(abs(X - 2) < 1e-8)

# Extract y samples
y_strong <- ys[, c(X_0_index, X_02_index)]
y_weak   <- ys[, c(X_0_index, X_2_index)]


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Strong
contour(
  surf_strong$x1, surf_strong$x2, surf_strong$pdf,
  nlevels = 25, col = heat.colors(25),
  xlab = expression(y[1]), ylab = expression(y[2]),
  main = paste0("Strong correlation\nk(0, 0.2) = ", round(Σ_strong[1,2], 2)),
  xlim = c(-2.7, 2.7), ylim = c(-2.7, 2.7)
)

points(
  y_strong, col = "black", pch = 16
)

# Weak
contour(
  surf_weak$x1, surf_weak$x2, surf_weak$pdf,
  nlevels = 25, col = heat.colors(25),
  xlab = expression(y[1]), ylab = expression(y[2]),
  main = paste0("Weak correlation\nk(0, 2) = ", round(Σ_weak[1,2], 2)),
  xlim = c(-2.7, 2.7), ylim = c(-2.7, 2.7)
)

points(
  y_weak, col = "black", pch = 16
)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP 4
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

GP_noise <- function(X1, y1, X2, kernel_func, sigma_noise) {
  n1 <- nrow(X1)

  # Covariance of noisy training observations
  Sigma11 <- kernel_func(X1, X1) + diag(sigma_noise^2, n1)

  # Cross-covariance
  Sigma12 <- kernel_func(X1, X2)

  # Solve Sigma11 * x = Sigma12 for x
  solved <- t(solve(Sigma11, Sigma12))

  # Posterior mean
  mu2 <- solved %*% y1

  # Posterior covariance
  Sigma22 <- kernel_func(X2, X2)
  Sigma2 <- Sigma22 - (solved %*% Sigma12)

  list(mean = mu2, covariance = Sigma2)
}

X1 <- matrix(c(-1, 0, 1), ncol = 1)
y1 <- c(0.5, -0.2, 0.3)
X2 <- matrix(seq(-2, 2, length.out = 50), ncol = 1)
sigma_noise <- 0.1
result <- GP_noise(X1, y1, X2, exponentiated_quadratic, sigma_noise)
mu2 <- result$mean
Sigma2 <- result$covariance

set.seed(123)  # for reproducibility
n1 <- nrow(X1)
sigma_noise <- 1
y1_noisy <- y1 + sigma_noise^2 * rnorm(n1)

result <- GP_noise(X1, y1_noisy, X2, exponentiated_quadratic, sigma_noise)
mu2 <- result$mean
Sigma2 <- result$covariance

sigma2 <- sqrt(diag(Sigma2))

ny <- 5  # number of posterior samples
set.seed(456)
L <- chol(Sigma2 + diag(1e-8, nrow(Sigma2)))  # add small jitter for stability
y2_samples <- matrix(NA, nrow = ny, ncol = nrow(X2))
for (i in 1:ny) {
  z <- rnorm(nrow(X2))
  y2_samples[i, ] <- mu2 + t(L) %*% z
}

f_sin <- function(x) {
  sin(x)
}

domain <- c(-6, 6)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# --- First plot: Posterior mean and uncertainty
plot(
  X2, mu2, type = "l", lwd = 2, col = "red",
  ylim = c(-3, 3), xlim = domain,
  xlab = "x", ylab = "y",
  main = "Distribution of posterior and prior data"
)

# Underlying true function
lines(X2, f_sin(X2), lty = 2, col = "blue")

# Uncertainty region
polygon(
  c(X2, rev(X2)),
  c(mu2 - 2 * sigma2, rev(mu2 + 2 * sigma2)),
  col = adjustcolor("red", alpha.f = 0.15), border = NA
)

# Training data
points(X1, y1_noisy, pch = 16)

legend(
  "topright",
  legend = c("True sin(x)", "Posterior mean", "Training points", "95% CI"),
  col = c("blue", "red", "black", adjustcolor("red", alpha.f = 0.15)),
  lty = c(2, 1, NA, NA),
  pch = c(NA, NA, 16, 15),
  pt.bg = c(NA, NA, "black", adjustcolor("red", alpha.f = 0.15)),
  bty = "n"
)

# --- Second plot: Posterior samples
matplot(
  X2, t(y2_samples),
  type = "l", lty = 1,
  col = 1:ny,
  xlab = "x", ylab = "y",
  main = "5 different function realizations from posterior",
  xlim = c(-6, 6), ylim = c(-3, 3)
)

points(X1, y1_noisy, pch = 16)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# STEP 4
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////


install.packages("DiceKriging")
library(DiceKriging)

set.seed(123)

# True function
f_sin <- function(x) sin(x)
X_domain <- c(-6, 6)

# Training data
X1 <- matrix(c(-5, -3, -1, 0, 1, 3, 5), ncol = 1)
sigma_noise <- 1
y1_true <- f_sin(X1)
y1_noisy <- y1_true + sigma_noise^2 * rnorm(length(X1))
y1_noisy

# Fit GP model
fit <- km(
  formula = ~X,
  design = data.frame(X = X1),
  response = y1_noisy,
  covtype = "gauss",
  nugget.estim = TRUE
)

# Predictions
X2 <- seq(-6, 6, length.out = 200)
pred <- predict(fit, newdata = data.frame(X = X2), type = "UK")

# Plot posterior mean and uncertainty
plot(
  X2, pred$mean, type = "l", lwd = 2, col = "red",
  ylim = c(-3, 3),
  xlab = "x", ylab = "y",
  main = "GP posterior with DiceKriging"
)
lines(X2, f_sin(X2), col = "blue", lty = 2)
polygon(
  c(X2, rev(X2)),
  c(pred$mean - 2 * pred$sd, rev(pred$mean + 2 * pred$sd)),
  col = adjustcolor("red", alpha.f = 0.15), border = NA
)
points(X1, y1_noisy, pch = 16)
legend(
  "topright",
  legend = c("True sin(x)", "Posterior mean", "Training points", "95% CI"),
  col = c("blue", "red", "black", adjustcolor("red", alpha.f = 0.15)),
  lty = c(2, 1, NA, NA),
  pch = c(NA, NA, 16, 15),
  pt.bg = c(NA, NA, "black", adjustcolor("red", alpha.f = 0.15)),
  bty = "n"
)

# Draw samples
set.seed(456)
ny <- 5
samples <- simulate(
  fit,
  nsim = ny,
  newdata = data.frame(X = X2),
  cond = TRUE
)

# Plot samples
matplot(
  X2, samples,
  type = "l", lty = 1,
  col = 1:ny,
  xlab = "x", ylab = "y",
  main = "5 different function realizations from posterior"
)








