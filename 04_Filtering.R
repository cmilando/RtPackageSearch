# Filtering

# INSTALL:
# remotes::install_github("dajmcdon/rtestim")

library(rtestim)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# GET GAMMA SERIAL INTERVAL
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# RtESTIM assumes that the serial interval is gamma
# here is code that will find a gamma for your discrete distribution

library(nloptr)

# serial interval, remove leading 0
x <- serial_interval_pmf[2:length(serial_interval_pmf)]

#
x[x == 0] <- 1e-6

#
my_eval <- function(params, data) {
  si_shape <- params[1]
  si_rate <- params[2]

  if (si_shape <= 0 || si_rate <= 0) return(1e10)  # penalize invalid params

  w <- sapply(1:(length(data)), function(x){
    pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
  })

  cor(w, data) * -1
}


# Initial parameter guesses
init_params <- c(shape = 2, rate = 10)

# Bounds (positive shape and rate)
lower_bounds <- c(1e-6, 1e-6)
upper_bounds <- c(100, 1000)

# Optimize using nloptr
fit <- nloptr(
  x0 = init_params,
  eval_f = my_eval,
  lb = lower_bounds,
  ub = upper_bounds,
  opts = list("algorithm" = "NLOPT_LN_SBPLX",
              "xtol_rel" = 1e-8,
              "maxeval" = 1000),
  data = x
)

# View estimated parameters
si_shape = fit$solution[1]
si_rate  = fit$solution[2]

w <- sapply(1:(length(serial_interval_pmf)), function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
})

# looks pretty good
plot(w)
lines(x, col = 'red')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

get_rtestim_output <- function(k) {

  # rt estimation
  rtestim <- cv_estimate_rt(
    dist_gamma      = c(si_shape, si_rate),
    observed_counts = reports_df$N,
    x               = reports_df$date,
    korder          = k,
    nsol            = 1000,
    maxiter         = 1e8
  )

  #approximate confidence bands
  rtestim_cb <- confband(rtestim, lambda = "lambda.1se")
  #lambda: the selected lambda. May be a scalar value,
  # or in the case of cv_poisson_rt objects, "lambda.min" or "lambda.max"

  # create dataframe
  # if you want to, shift backwards by the seeding_time, aka the sum
  # of the means of the delay distributions
  plot_rtestim <- data.frame(
    model = paste0("k = ",k),
    ## *****
    ## date = reports_df$date - seeding_time,
    date = reports_df$date,
    ###
    Rt = rtestim_cb$fit,
    Rt_lb = rtestim_cb$`2.5%`,
    Rt_ub = rtestim_cb$`97.5%`
  )
  return(plot_rtestim)
}

df1 <- get_rtestim_output(1)
df3 <- get_rtestim_output(5)
rtestim_df <- rbind(df1, df3)

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

plot_rt1 <- ggplot(rtestim_df) +
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
  scale_color_discrete(name = 'Filter degree') +
  scale_fill_discrete(name = 'Filter degree') +
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

ggsave('img/Filtering.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/Filtering.RDS")

