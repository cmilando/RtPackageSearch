# Fixed Sliding Windows

# install.packages('EpiEstim',
#                  repos = c('https://mrc-ide.r-universe.dev',
#                            'https://cloud.r-project.org'))
library(EpiEstim)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

setnames(reports_df, "N", "I")

get_epistim_output <- function(window_size) {

  # THEN SET SLIDING WINDOW SIZE
  estimation_window <- window_size

  t_start <- seq(2, nrow(reports_df) - estimation_window)
  t_end <- t_start + estimation_window

  reports_df$date_int <- 1:nrow(reports_df)

  # THEN ESTIMATE R(t) USING SLIDING WINDOWS
  getR <- EpiEstim::estimate_R(
    incid = reports_df,
    method = "non_parametric_si",
    config = make_config(list(
      si_distr = serial_interval_pmf,
      t_start = t_start,
      t_end = t_end
    )),
    backimputation_window = 15
  )

  # **********
  # INCLUDE THE DECONVOLUTION
  # getR$R$date_int <- getR$R$t_end - seeding_time
  getR$R$date_int <- getR$R$t_end
  # **********

  # **********
  # MOVE TO WINDOW CENTER
  getR$R$date_int <- getR$R$t_end - floor(estimation_window/2)
  # **********

  #
  EpiEstim_R <- getR$R[, c('t_start', 't_end', 'date_int', 'Median(R)',
                           'Quantile.0.05(R)', 'Quantile.0.95(R)')]

  names(EpiEstim_R) <- c('t_start', 't_end', "date_int", "Rt", "Rt_lb", "Rt_ub")

  EpiEstim_R$model <- paste0(window_size, " days")
  EpiEstim_R$date_int <- as.integer(EpiEstim_R$date_int)
  reports_df <- data.frame(reports_df)

  EpiEstim_R <- merge(EpiEstim_R, reports_df[, c('date', 'date_int')])

  return(EpiEstim_R)

}


EpiEstim_R1 <- get_epistim_output(2)
EpiEstim_R2 <- get_epistim_output(9)

EpiEstim_R_full <- rbind(EpiEstim_R1, EpiEstim_R2)

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

plot_rt1 <- ggplot(EpiEstim_R_full) +
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
  scale_color_discrete(name = 'Sliding window') +
  scale_fill_discrete(name = 'Sliding window') +
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
ggsave('img/FixedSlidingWindow.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/FixedSlidingWindow.RDS")
