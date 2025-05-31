# Random Walk


# install.packages("EpiNow2",
#                  repos = c("https://epiforecasts.r-universe.dev",
#                            getOption("repos")))

#
library(EpiNow2)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

##
gi_pmf         <- NonParametric(pmf = generation_interval_pmf)
infect_to_test <- NonParametric(pmf = infect_to_test_pmf)
sym_report_delay_pmf <- NonParametric(pmf = reporting_delay_pmf)

mean(gi_pmf) + mean(infect_to_test) + mean(sym_report_delay_pmf)

setnames(reports_df, 'N', 'confirm')

## -------------------
# Including all the options so I can see how  this works
# see defaults here: https://epiforecasts.io/EpiNow2/reference/obs_opts.html
get_EpiNow2output <- function(rw) {

  res_epinow <- epinow(
    data                 = reports_df,
    generation_time      = generation_time_opts(gi_pmf),
    delays               = delay_opts(infect_to_test),
    truncation           = trunc_opts(sym_report_delay_pmf),
    rt                   = rt_opts(rw = rw),
    gp                   = NULL, # turn off the GP model here!
    backcalc             = backcalc_opts(prior = 'reports'),
    stan                 = stan_opts(chains = 4, cores = 4),
    obs                  = obs_opts(), # ok to use these defaults
    forecast             = forecast_opts(),
    CrIs                 = c(0.2, 0.5, 0.9)
  )

  R_df <- subset(res_epinow$estimates$summarised, variable == 'R')

  R_df$model <- paste0(rw, " days")
  R_df$Rt <- R_df$median
  R_df$Rt_lb <- R_df$lower_90
  R_df$Rt_ub <- R_df$upper_90
  return(R_df)
}

R_df <- get_EpiNow2output(2)
R_df2 <- get_EpiNow2output(6)

##
R_df_all <- rbind(R_df, R_df2)
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

plot_rt1 <- ggplot(R_df_all) +
  ##
  theme_classic2() +
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
  scale_color_discrete(name = 'Random walk') +
  scale_fill_discrete(name = 'Random walk') +
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
  annotate('text', x = first_day + 35,
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

# dev.size()
ggsave('img/RandomWalk.png', width = 6.5, height = 2)
saveRDS(plot_rt1, "img/RandomWalk.RDS")

