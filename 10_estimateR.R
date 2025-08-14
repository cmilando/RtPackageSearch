
# https://github.com/covid-19-Re/estimateR

# library(devtools)
# install_github("covid-19-Re/estimateR")
library(estimateR)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////


# unclear if there is functionality at present to use non-parametric distributions
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
onset_to_report <- list(name="gamma",
                        shape = shape_onset_to_report,
                        scale = scale_onset_to_report)

shape_incubation = 3.2
scale_incubation = 1.3
incubation <- list(name="gamma",
                   shape = shape_incubation,
                   scale = scale_incubation)


mean_serial_interval = 4.8
std_serial_interval = 2.3

estimation_window = 3

N_bootstrap_replicates <- 100

estimates <- get_block_bootstrapped_estimate(
  reports_df$N,
  N_bootstrap_replicates = N_bootstrap_replicates,
  smoothing_method = "LOESS",
  deconvolution_method = "Richardson-Lucy delay distribution",
  estimation_method = "EpiEstim sliding window",
  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
  combine_bootstrap_and_estimation_uncertainties = TRUE,
  delay = list(incubation, onset_to_report),
  estimation_window = estimation_window,
  mean_serial_interval = mean_serial_interval,
  std_serial_interval  = std_serial_interval,
  output_Re_only = FALSE,
  ref_date = reports_df$date[1],
  time_step = "day"
)


#
estimateR_R <- estimates
estimateR_R$Rt <- estimateR_R$Re_estimate
estimateR_R$Rt_lb <- estimateR_R$Re_lowHPD
estimateR_R$Rt_ub <- estimateR_R$Re_highHPD

estimateR_R$model <- 'Standard'

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

plot_rt1 <- ggplot(estimateR_R) +
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
  scale_color_discrete(name = 'estimateR') +
  scale_fill_discrete(name = 'estimateR') +
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
ggsave('img/estimateR.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/estimateR.RDS")
