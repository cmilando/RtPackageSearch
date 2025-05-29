# Fixed Sliding Windows

# install.packages('EpiEstim',
#                  repos = c('https://mrc-ide.r-universe.dev',
#                            'https://cloud.r-project.org'))
library(EpiEstim)

#
library(ggplot2)
library(ggpubr)
library(patchwork)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# AGGREGATED REPORTED CASE INCIDENCE DATA
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# create some Aggregated incidence data
create_agg_data <- function() {
  set.seed(123)
  days      <- 1:80
  amplitude <- 0.5          # Maximum height of the sine wave
  period    <- 40            # Total time period
  frequency <- 2          # Two peaks in the 50-day period
  sine_wave <- amplitude * sin(2 * pi * frequency * days / period)
  seq1      <- exp(sine_wave[15:40]) * 1000
  init      <- exp(runif(15)) * 100
  lambda    <- c(init[1:8], seq1[1:20], init[9:15], seq1[21:26]*0.65)
  cases     <- sapply(lambda, function(l) rpois(1, l + exp(rnorm(1, 0, 2))))
  return(cases)
}

incidence_df <- EpiNow2::example_confirmed
incidence_df <- incidence_df[c(10:50), ]
incidence_df$confirm <- create_agg_data()
incidence_df$I <- incidence_df$confirm

day1 <- incidence_df$date[1]
plot(incidence_df$confirm)
# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# GENERATION TIME and SERIAL INTERVAL
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

## Has to start with 0 !
generation_interval_pmf <- c(
  0, 0.0610, 0.1540, 0.2198, 0.2178, 0.1536, 0.1122, 0.0486, 0.0224,
  0.0078, 0.0022, 0.0004, 0.0002
)

## ASSUME THAT SERIAL INTERVAL = GENERATION INTERVAL
serial_interval_pmf <- generation_interval_pmf

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# DELAYS
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# 1) Infection to taking a test
# has to start with 0 !
infect_to_test_pmf <- c(
  0, 0.0001, 0.0030, 0.1422, 0.2714, 0.2664,
  0.1832, 0.0846, 0.0340, 0.0136, 0.0014, 0.0001
)

# 2) time from taking a test to it getting into a state database
reporting_delay_pmf <- c(
  0.3786, 0.3724, 0.1662, 0.0622, 0.0166, 0.0034, 0.0006
)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# ADD RIGHT TRUNCATION TO DATA BASED ON REPORTING DELAY
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

tail_Px <- rev(cumsum(reporting_delay_pmf))
px_len <- length(tail_Px)

incidence_df$Px = 1
cases_len <- nrow(incidence_df)

incidence_df$Px[(cases_len - px_len + 1):cases_len] <- tail_Px

incidence_df$confirm_true <- incidence_df$confirm
incidence_df$confirm      <- incidence_df$confirm * incidence_df$Px
incidence_df$confirm      <- as.integer(incidence_df$confirm)
incidence_df$date_int     <- 1:nrow(incidence_df)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# EPIESTIM DOES NOT DO NOW-CASTING
# but here is where you would do it if it did

# FIRST DECONVOLVE TO INFECTIONS FROM REPORTS
set.seed(123)

joint_delay <- sapply(1:1e4, function(i) {
  sample(0:(length(infect_to_test_pmf) - 1), size = 1,
         prob = infect_to_test_pmf) +
  sample(0:(length(reporting_delay_pmf) - 1), size = 1,
         prob = reporting_delay_pmf)
})

D <- data.frame(i = 0:max(joint_delay), n = 0, Px = NA)
for(i in 0:max(joint_delay)) {
  D$n[i+1] = length(which(joint_delay == i))
}

D$Px = D$n / sum(D$n)
D_delay <- round(weighted.mean(x = 1:length(D$Px), w = D$Px))

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

plot_report <- ggplot(subset(incidence_df, date >= day1)) +
  theme_classic2() +
  geom_vline(xintercept = c(as.Date("2020-03-30") + 0.5,
                            as.Date("2020-04-11") + 0.5),
             linetype = '41') +
  geom_vline(xintercept = as.Date("2020-04-11") -
               length(reporting_delay_pmf) + 0.5,
             linetype = '21', color = 'pink') +
  geom_col(aes(x = date, y = confirm),
           fill = grey(0.75), width = 0.5) +
  coord_cartesian(xlim = c(as.Date("2020-03-02"),
                           as.Date("2020-04-18")),
                  ylim = c(0, report_ymax)) +
  ylab("Reported cases") +
  xlab('Date') +
  annotate('text', x = as.Date("2020-03-05"),
           y = report_ymax, label = 'Historical period', size = 3) +
  annotate('text', x = as.Date("2020-04-02"),
           y = report_ymax, label = 'Nowcast', size = 3) +
  annotate('text', x = as.Date("2020-04-05"),
           y = report_ymax, label = 'Right-\ntruncation', size = 2,
           vjust = 0.80, hjust = 0, color = scales::muted('red'))+
  annotate('text', x = as.Date("2020-04-14"),
           y = report_ymax, label = 'Forecast', size = 3)


plot_rt1 / (plot_report) + plot_layout(guides = "collect")

# dev.size()
ggsave('img/FixedSlidingWindow.png', width = 6*1.5, height = 3.5*1.3)
