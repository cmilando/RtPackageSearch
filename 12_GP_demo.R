# implemented Guassian Process model
# Adapted from these sources:
# https://peterroelants.github.io/posts/gaussian-process-tutorial/
# https://rpubs.com/ncahill_stat/889056
# https://www.rdocumentation.org/packages/brms/versions/1.10.2/topics/gp
# https://medium.com/@leungcheuk209/gaussian-process-with-r-examples-440d273629a9
# https://www.jstatsoft.org/article/view/v109i05

library(tidyverse)
library(patchwork)
library(ggpubr)

library(greta)
library(greta.gp)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# SIMULATE DATA
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

set.seed(123)

X1 <- matrix(seq(from = 0, to = 20, by = 1), ncol = 1)
X2 <- matrix(seq(from = min(X1), to = max(X1), length.out = 1000),
             ncol = 1)

y1 <- matrix(sin(X1) + rnorm(length(X1), 0, 0.5), ncol= 1)

ggplot() + theme_classic2() +
  geom_point(aes(x = X1, y = y1), size = 1) +
  geom_line(aes(x = X2, y = sin(X2)),
            color = 'blue', linetype = 'dashed')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# Exponentiated quadratic kernel
exponentiated_quadratic_kernel <- function(xa, xb,
                                           amplitude = 0.5,
                                           lengthscale = 0.5) {
  # Ensure column vectors
  xa <- as.matrix(xa)
  xb <- as.matrix(xb)

  # Compute squared distances
  dists <- as.matrix(dist(rbind(xa, xb)))^2
  n_a <- nrow(xa)
  n_b <- nrow(xb)
  sq_dists <- dists[1:n_a, (n_a + 1):(n_a + n_b)]

  # Compute kernel
  k <- amplitude^2 * exp(-0.5 * sq_dists / lengthscale^2)
  return(k)
}


GP <- function(X1, y1, X2, kernel_func) {

  # Calculate the posterior mean and covariance matrix for y2
  # based on the corresponding input X2, the observations (y1, X1),
  # and the prior kernel function.

  # Kernel of the observations
  Sigma11 <- kernel_func(X1, X1)

  # Kernel of observations (x1) vs to-predict (x2)
  Sigma12 <- kernel_func(X1, X2)

  # Solve
  solved <- t(solve(Sigma11, Sigma12))

  # Compute Posterior mean
  mu2 <- solved %*% y1

  # Compute the posterior covariance
  Sigma22 <- kernel_func(X2, X2)
  Sigma2 <- Sigma22 - (solved %*% Sigma12)

  # return mean and covariance
  list(mean = mu2, covariance = Sigma2)
}

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# CALCULATE
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

result <- GP(X1, y1, X2, exponentiated_quadratic_kernel)

mu2 <- result$mean

Sigma2 <- result$covariance
sigma2 <- sqrt(ifelse(diag(Sigma2)<0, 0, diag(Sigma2)))

get_y2 <- function(mu2, Sigma2, ny = 5) {
  ny <- 5  # number of posterior samples
  L <- chol(Sigma2 + diag(1e-8, nrow(Sigma2)))  # add small jitter for stability
  y2_samples <- matrix(NA, nrow = ny, ncol = nrow(X2))
  for (i in 1:ny) {
    z <- rnorm(nrow(X2))
    y2_samples[i, ] <- mu2 + t(L) %*% z
  }
  return(y2_samples)
}

y2_samples <- get_y2(mu2 = mu2, Sigma2 = Sigma2)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

p0 <- ggplot() + theme_classic2() +
  # geom_ribbon(aes(x = X2,
  #                 ymin = mu2 - 2 * sigma2,
  #                 ymax = mu2 + 2 * sigma2),
  #             fill = adjustcolor("red", alpha.f = 0.15)) +
  # geom_line(aes(x = X2, y = mu2), color = 'red') +
  geom_line(aes(x = X2, y = sin(X2)), color = 'blue',
            linetype = 'dashed') +
  geom_point(aes(x = X1, y = y1), size = 1.5) + xlab("X") + ylab("Y") +
  ggtitle("a.")

# --- First plot: Posterior mean and uncertainty
p1 <- ggplot() + theme_classic2() +
  geom_ribbon(aes(x = X2,
                  ymin = mu2 - 2 * sigma2,
                  ymax = mu2 + 2 * sigma2),
              fill = adjustcolor("red", alpha.f = 0.15)) +
  geom_line(aes(x = X2, y = mu2), color = 'red') +
  # geom_line(aes(x = X2, y = sin(X2)), color = 'blue',
  #           linetype = 'dashed') +
  geom_point(aes(x = X1, y = y1), size = 1.5) + xlab("X") + ylab("Y") +
  ggtitle("c.")

# --- Second plot: Posterior samples
p2 <- ggplot() + theme_classic2() +
  geom_line(aes(x = X2, y = y2_samples[1, ]), color = 'orange') +
  geom_line(aes(x = X2, y = y2_samples[2, ]), color = 'orange') +
  geom_line(aes(x = X2, y = y2_samples[3, ]), color = 'orange') +
  geom_line(aes(x = X2, y = y2_samples[4, ]), color = 'orange') +
  geom_line(aes(x = X2, y = y2_samples[5, ]), color = 'orange') +
  geom_point(aes(x = X1, y = y1), size = 1.5) +
  xlab("X") + ylab("Y") +  ggtitle("b.")

p0 / p2 / p1
ggsave("gp_supp_fig1.png", width = 6, height = 6)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# EXTEND TO NOISY DATA
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# hyperparameters
rbf_var <- lognormal(0, 1)  #
rbf_len <- lognormal(0, 1)  #
obs_sd <- lognormal(0, 1)   #

# kernel & GP
kernel <- rbf(rbf_len, rbf_var) + bias(1)
f <- gp(X1, kernel)

# likelihood
distribution(y1) <- normal(f, obs_sd)

# prediction
f_plot <- project(f, X2)

# fit the model by Hamiltonian Monte Carlo
m <- model(f_plot)
draws <- mcmc(m, chains = 4, n_samples = 8000)

dim(draws[[1]])

greta_med = unname(apply(draws[[1]], 2, function(x) quantile(x, probs = 0.5)))
greta_lb  = unname(apply(draws[[1]], 2, function(x) quantile(x, probs = 0.025)))
greta_ub  = unname(apply(draws[[1]], 2, function(x) quantile(x, probs = 0.975)))

p3 <- ggplot() + theme_classic2() +
  geom_ribbon(aes(x = X2, ymin = greta_lb, ymax = greta_ub),
              fill = 'lavender', alpha = 0.5) +
  geom_line(aes(x = X2, y = greta_med), color ='purple') +
  geom_line(aes(x = X2, y = sin(X2)), color = 'blue',
            linetype = 'dashed') +
  geom_point(aes(x = X1, y = y1)) +
  xlab("X") + ylab("Y")

p3

ggsave("gp_supp_fig2.png", width = 6, height = 2)
