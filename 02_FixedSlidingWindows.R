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

# EPIESTIM DOES NOT DO NOW-CASTING
# but here is where you would do it if it did

# THEN SET SLIDING WINDOW SIZE
estimation_window <- 2

t_start <- seq(2, nrow(incidence_df) - estimation_window)
t_end <- t_start + estimation_window

# THEN ESTIMATE R(t) USING SLIDING WINDOWS
getR <- EpiEstim::estimate_R(
  incid = incidence_df,
  method = "non_parametric_si",
  config = make_config(list(
    si_distr = serial_interval_pmf,
    t_start = t_start,
    t_end = t_end
  )),
  backimputation_window = 5
)

# INCLUDE THE DECONVOLUTION
getR$R$date_int <- getR$R$t_end - D_delay

#
EpiEstim_R <- getR$R[, c('t_start', 't_end', 'date_int', 'Median(R)',
                         'Quantile.0.05(R)', 'Quantile.0.95(R)')]

names(EpiEstim_R) <- c('t_start', 't_end', "date_int", "Rt", "Rt_lb", "Rt_ub")

EpiEstim_R$model <- '2 days'
EpiEstim_R$date_int <- as.integer(EpiEstim_R$date_int)
incidence_df <- data.frame(incidence_df)

EpiEstim_R <- merge(EpiEstim_R, incidence_df[, c('date','date_int')])
EpiEstim_R_full <- EpiEstim_R

# Repeat again to get another estimate
estimation_window <- 9
t_start <- seq(2, nrow(incidence_df) - estimation_window)
t_end <- t_start + estimation_window

getR <- EpiEstim::estimate_R(
  incid = incidence_df ,
  method = "non_parametric_si",
  config = make_config(list(
    si_distr = serial_interval_pmf,
    t_start = t_start,
    t_end = t_end
  )),
  backimputation_window = 5
)

# INCLUDE DECONVOLUTION
getR$R$date_int <- getR$R$t_end - D_delay

#
EpiEstim_R <- getR$R[, c('t_start', 't_end', 'date_int', 'Median(R)',
                         'Quantile.0.05(R)', 'Quantile.0.95(R)')]

names(EpiEstim_R) <- c('t_start', 't_end', "date_int", "Rt", "Rt_lb", "Rt_ub")

EpiEstim_R$model <- '9 days'
EpiEstim_R$date_int <- as.integer(EpiEstim_R$date_int)
EpiEstim_R <- merge(EpiEstim_R, incidence_df[, c('date','date_int')])

EpiEstim_R_full <- rbind(EpiEstim_R, EpiEstim_R_full)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

report_ymax <- 4000
rt_ymax <- 4.0

plot_rt1 <- ggplot(subset(EpiEstim_R_full, date >= day1)) +
  theme_classic2() +
  geom_hline(yintercept = 1, linetype = '11') +
  geom_vline(xintercept = c(as.Date("2020-03-30") + 0.5,
                            as.Date("2020-04-11") + 0.5),
             linetype = '41') +
  geom_vline(xintercept = as.Date("2020-04-11") -
               length(reporting_delay_pmf) + 0.5,
             linetype = '21', color = 'pink') +
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model), linewidth = 0.25) +
  geom_point(aes(x = date, y = Rt, color = model), size = 0.5) +
  scale_color_discrete(name = 'Sliding window size') +
  scale_fill_discrete(name = 'Sliding window size') +
  coord_cartesian(xlim = c(as.Date("2020-03-02"),
                           as.Date("2020-04-18")),
                  ylim = c(0, rt_ymax)) +
  ylab(expression(R[t])) +
  annotate('text', x = as.Date("2020-03-05"),
           y = rt_ymax, label = 'Historical period', size = 3) +
  annotate('text', x = as.Date("2020-04-02"),
           y = rt_ymax, label = 'Nowcast', size = 3) +
  annotate('text', x = as.Date("2020-04-05"),
           y = rt_ymax, label = 'Right-\ntruncation', size = 2,
           vjust = 0.80, hjust = 0, color = scales::muted('red'))+
  annotate('text', x = as.Date("2020-04-14"),
           y = rt_ymax, label = 'Forecast', size = 3) +
  xlab(NULL)

plot_rt1

# plot_report <- ggplot(subset(incidence_df, date >= day1)) +
#   theme_classic2() +
#   geom_vline(xintercept = c(as.Date("2020-03-30") + 0.5,
#                             as.Date("2020-04-11") + 0.5),
#              linetype = '41') +
#   geom_vline(xintercept = as.Date("2020-04-11") -
#                length(reporting_delay_pmf) + 0.5,
#              linetype = '21', color = 'pink') +
#   geom_col(aes(x = date, y = confirm),
#            fill = grey(0.75), width = 0.5) +
#   coord_cartesian(xlim = c(as.Date("2020-03-02"),
#                            as.Date("2020-04-18")),
#                   ylim = c(0, report_ymax)) +
#   ylab("Reported cases") +
#   xlab('Date') +
#   annotate('text', x = as.Date("2020-03-05"),
#            y = report_ymax, label = 'Historical period', size = 3) +
#   annotate('text', x = as.Date("2020-04-02"),
#            y = report_ymax, label = 'Nowcast', size = 3) +
#   annotate('text', x = as.Date("2020-04-05"),
#            y = report_ymax, label = 'Right-\ntruncation', size = 2,
#            vjust = 0.80, hjust = 0, color = scales::muted('red'))+
#   annotate('text', x = as.Date("2020-04-14"),
#            y = report_ymax, label = 'Forecast', size = 3)
#
#
# plot_rt1 / (plot_report) + plot_layout(guides = "collect")
#
# # dev.size()
# ggsave('img/FixedSlidingWindow.png', width = 6*1.5, height = 3.5*1.3)
