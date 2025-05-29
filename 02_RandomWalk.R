# Random Walk

#
library(EpiNow2)

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

## Has to start with 0 !
gi_pmf <- NonParametric(pmf = generation_interval_pmf)

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

infect_to_test <- NonParametric(pmf = infect_to_test_pmf)

# 2) time from taking a test to it getting into a state database
reporting_delay_pmf <- c(
  0.3786, 0.3724, 0.1662, 0.0622, 0.0166, 0.0034, 0.0006
)

sym_report_delay_pmf <- NonParametric(pmf = reporting_delay_pmf)

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

## -------------------
rw = 2

res_epinow <- epinow(
  data            = incidence_df[, c('date', 'confirm')],
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  rt              = rt_opts(rw = rw),
  truncation      = trunc_opts(sym_report_delay_pmf),
  backcalc        = backcalc_opts(prior = 'reports'),
  stan            = stan_opts(chains = 4, cores = 4)
)

plot(res_epinow)

R_df <- subset(res_epinow$estimates$summarised, variable == 'R')
R_df$date_int <- 1:nrow(R_df) - 1
R_df$model <- paste0(rw, " days")
R_df$Rt <- R_df$median
R_df$Rt_lb <- R_df$lower_90
R_df$Rt_ub <- R_df$upper_90

cases_df <- subset(res_epinow$estimates$summarised, variable == 'reported_cases')
cases_df$date_int <- 1:nrow(R_df) - 1
cases_df$model <- paste0(rw, " day")
cases_df$reports <- cases_df$median
cases_df$reports_lb <- cases_df$lower_90
cases_df$reports_ub <- cases_df$upper_90

## -------------------
rw = 7

res_epinow2 <- epinow(
  data            = incidence_df[, c('date', 'confirm')],
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  rt              = rt_opts(rw = rw),
  truncation      = trunc_opts(sym_report_delay_pmf),
  backcalc        = backcalc_opts(prior = 'reports'),
  stan            = stan_opts(chains = 4, cores = 4)
)

R_df2 <- subset(res_epinow2$estimates$summarised, variable == 'R')
R_df2$date_int <- 1:nrow(R_df2) - 1
R_df2$model <- paste0(rw, " days")
R_df2$Rt <- R_df2$median
R_df2$Rt_lb <- R_df2$lower_90
R_df2$Rt_ub <- R_df2$upper_90

cases_df2 <- subset(res_epinow2$estimates$summarised, variable == 'reported_cases')
cases_df2$date_int <- 1:nrow(R_df2) - 1
cases_df2$model <- paste0(rw, " days")
cases_df2$reports <- cases_df2$median
cases_df2$reports_lb <- cases_df2$lower_90
cases_df2$reports_ub <- cases_df2$upper_90

##
R_df_all <- rbind(R_df, R_df2)
cases_df_all <- rbind(cases_df, cases_df2)
lastday <- max(R_df_all$date)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

report_ymax <- 4000
rt_ymax <- 4.0

plot_rt1 <- ggplot(subset(R_df_all, date >= day1)) +
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
  scale_color_discrete(name = 'Random walk size') +
  scale_fill_discrete(name = 'Random walk size') +
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
  geom_ribbon(data = cases_df_all,
              aes(x = date,
                  ymin = reports_lb,
                  ymax = reports_ub,fill = model),
              alpha = 0.25, show.legend = F) +
  geom_line(data = cases_df_all,
            aes(x = date, y = reports, color = model),
            show.legend = F)+
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
ggsave('img/RandomWalk.png', width = 6*1.5, height = 3.5*1.3)
