# B-Splines

library(EpiLPS)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

get_epilps_output <- function(k) {

  # discrete dist
  si_spec <- Idist(probs = serial_interval_pmf)

  #
  LPSfit <- estimR(incidence = reports_df$N,
                   si = si_spec$pvec,
                   K = k)

  # create dataframe
  # shift backwards by the seeding_time, aka the sum
  # of the means of the delay distributions
  plot_epilps <- data.frame(
    model = paste0("k = ",k),
    date = reports_df$date - seeding_time,
    Rt = LPSfit$RLPS$Rq0.50,
    Rt_lb = LPSfit$RLPS$Rq0.025,
    Rt_ub = LPSfit$RLPS$Rq0.975
  )
  return(plot_epilps)
}

df1 <- get_epilps_output(30)
df3 <- get_epilps_output(9)
epilps_df <- rbind(df1, df3)

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

plot_rt1 <- ggplot(epilps_df) +
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
  scale_color_discrete(name = '# of B-Splines') +
  scale_fill_discrete(name = '# of B-Splines') +
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

ggsave('img/BSplines.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/BSplines.RDS")

