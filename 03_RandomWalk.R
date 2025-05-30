# Random Walk

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
rw = 2

res_epinow <- epinow(
  data            = reports_df,
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  truncation      = trunc_opts(sym_report_delay_pmf),
  rt              = rt_opts(rw = rw),
  backcalc        = backcalc_opts(prior = 'reports'),
  stan            = stan_opts(chains = 4, cores = 4)
)

plot(res_epinow)

R_df <- subset(res_epinow$estimates$summarised, variable == 'R')
R_df$type
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
rw = 6

res_epinow2 <- epinow(
  data            = reports_df,
  generation_time = generation_time_opts(gi_pmf),
  delays          = delay_opts(infect_to_test),
  rt              = rt_opts(rw = rw),
  truncation      = trunc_opts(sym_report_delay_pmf),
  backcalc        = backcalc_opts(prior = 'reports'),
  stan            = stan_opts(chains = 4, cores = 4)
)

plot(res_epinow2)

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

plot_rt1 <- ggplot(subset(R_df_all)) +
  theme_classic2() +
  geom_hline(yintercept = 1, linetype = '11') +
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model), linewidth = 0.25) +
  geom_point(aes(x = date, y = Rt, color = model), size = 0.5) +
  scale_color_discrete(name = 'Random walk size') +
  scale_fill_discrete(name = 'Random walk size') +
  ylab(expression(R[t])) +
  xlab(NULL)

plot_rt1

# dev.size()
ggsave('img/RandomWalk.png', width = 6*1.5, height = 3.5*1.3)
