

# https://github.com/phac-nml-phrsd/ern

# https://cran.r-project.org/web/packages/ern/vignettes/est-rt.html

library(ern)
library(nloptr)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

names(reports_df)[2] <- 'value'

# NOTE in ERN we CANNOT DEFINE A non-parametric PMF, so will need
## to fit values for each one
get_gamma_params <- function(pmf_vec, seed=123) {
  set.seed(seed)
  shape_guess = 1
  rate_guess = 1
  x0 <- c(shape_guess, rate_guess)

  eval_f <- function(x0) {
    # get the test pmf
    sam = rgamma(1000, shape = x0[1], rate = x0[2])
    test_vec = table(cut(sam, breaks = 0:(length(pmf_vec)-1)))/1000
    while(length(test_vec) != length(pmf_vec)) {
      test_vec <- c(test_vec, 0)
    }
    # compare
    sum((pmf_vec - test_vec)^2)
  }

  ss <- nloptr(x0, eval_f, lb = c(0, 0), ub = c(100, 100),
               opts = list(algorithm = 'NLOPT_LN_COBYLA'))

  sam = rgamma(1000, shape = ss$solution[1], rate = ss$solution[2])
  test_vec = table(cut(sam, breaks = 0:(length(pmf_vec)-1)))/1000
  plot(x = 0:(length(pmf_vec)-1), y = c(test_vec, 0), type = 'l',
       ylim = c(0, max(pmf_vec, test_vec)))
  points(x = 0:(length(pmf_vec)-1), y = pmf_vec, col = 'red')

  return(ss$solution)
}


get_lognorm_params <- function(pmf_vec, seed = 123) {
  set.seed(seed)
  mean_guess = 1
  sd_guess = 1
  x0 <- c(mean_guess, sd_guess)

  eval_f <- function(x0) {
    # get the test pmf
    sam = exp(rnorm(1000, mean = x0[1], sd = x0[2]))
    test_vec = table(cut(sam, breaks = 0:(length(pmf_vec)-1)))/1000
    while(length(test_vec) != length(pmf_vec)) {
      test_vec <- c(test_vec, 0)
    }
    # compare
    sum((pmf_vec - test_vec)^2)
  }

  ss <- nloptr(x0, eval_f, lb = c(0, 0), ub = c(100, 100),
         opts = list(algorithm = 'NLOPT_LN_COBYLA'))

  sam = exp(rnorm(1000, mean = ss$solution[1], sd = ss$solution[2]))
  test_vec = table(cut(sam, breaks = 0:(length(pmf_vec)-1)))/1000
  plot(x = 0:(length(pmf_vec)-1), y = c(test_vec, 0), type = 'l',
       ylim = c(0, max(pmf_vec, test_vec)))
  points(x = 0:(length(pmf_vec)-1), y = pmf_vec, col = 'red')

  return(ss$solution)
}


#
# looks gamma so estimate
ss <- get_gamma_params(reporting_delay_pmf)
# the naming here doens't quite make sense given that this is gamma
dist.repdelay = ern::def_dist(
  dist = 'gamma',
  mean = ss[1],
  mean_sd = .1,
  sd = ss[2],
  sd_sd = .01,
  max = length(reporting_delay_pmf) - 1
)

# define reporting fraction -- assumed 100%
dist.repfrac = ern::def_dist(
  dist = "unif",
  min = 0.99,
  max = 1.0
)

# define incubation period
# this is some fraction of infect_to_test
# use the vignette version, unclear how to discretize this
dist.incub0 = ern::def_dist(
  dist     = "gamma",
  mean     = 3.49,
  mean_sd  = 0.1477,
  shape    = 8.5,
  shape_sd = 1.8945,
  max      = 11
)

# define generation interval
ss <- get_lognorm_params(generation_interval_pmf) # seems close enough
dist.gi = ern::def_dist(
  dist     = "lnorm",
  meanlog     = exp(ss[1]),
  meanlog_sd = 0.1,
  sdlog       = exp(ss[2]),
  sdlog_sd = 0.01,
  max      = length(generation_interval_pmf) - 1
)

r.estim = estimate_R_cl(
  cl.data      = reports_df,
  dist.repdelay = dist.repdelay,
  dist.repfrac  = dist.repfrac,
  dist.incub    = dist.incub0,
  dist.gi       = dist.gi
)

g = plot_diagnostic_cl(r.estim)
plot(g)

#
ern_R <- r.estim$R
ern_R$Rt <- ern_R$mean
ern_R$Rt_lb <- ern_R$lwr
ern_R$Rt_ub <- ern_R$upr

ern_R$model <- 'Standard'

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

rt_max <- 3
first_day <- min(reports_df$date)
last_day <- max(reports_df$date)
nowcast_start    = last_day - seeding_time
forecast_window  = last_day + 15

plot_rt1 <- ggplot(ern_R) +
  ##
  theme_classic2() +
  # ##
  annotate('rect',
           xmin = as.Date('2020-03-13') - 0.5,
           xmax = as.Date('2020-03-25') + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = 'lightyellow',
           alpha = 0.75,
           color = 'white') +
  annotate('rect',
           xmin = as.Date('2020-04-16') - 0.5,
           xmax = as.Date('2020-04-23') + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = 'lightyellow',
           alpha = 0.75,
           color = 'white') +
  #
  geom_hline(yintercept = 1, linetype = '11') +
  # ##
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model),
            linewidth = 0.25, show.legend = T) +
  geom_point(aes(x = date, y = Rt, color = model),
             size = 0.5,shape = 1,
             show.legend = T) +
  scale_color_discrete(name = 'ern') +
  scale_fill_discrete(name = 'ern') +
  ##
  coord_cartesian(xlim = c(first_day,
                           forecast_window),
                  ylim = c(0, rt_max+0.15), expand = F) +
  ylab(expression(R[t])) +
  xlab(NULL) +
  ##
  scale_x_date(breaks = '2 week',
               date_minor_breaks = "1 weeks",
               date_labels = "%b %d") +
  ##
  annotate('text', x = first_day + 55,
           y = rt_max, label = 'Historical period',
           size = 2) +
  annotate('text', x = nowcast_start + 6,
           y = rt_max, label = 'Nowcasting', size = 2) +
  annotate('text', x = last_day + 6,
           y = rt_max, label = 'Forecasting', size = 2) +
  ##
  geom_vline(xintercept = c(nowcast_start + 0.5,
                            last_day + 0.5),
             linewidth = 1.25,
             alpha = 0.5,
             linetype = 'solid', color = 'white') +
  geom_vline(xintercept = c(nowcast_start + 0.5),
             linetype = '22') +
  geom_vline(xintercept = c(
    last_day + 0.5),
    linetype = '41') +
  annotate('label',
           x = forecast_window - 14.5,
           size = 2,
           #label.size = 0,
           color = 'black',
           y = 2.2,
           angle = 90,
           label = 'PRESENT')

plot_rt1
ggsave('img/ern.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/ern.RDS")
