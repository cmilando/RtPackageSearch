# B-Splines

library(EpiLPS)

#
source('01_ReportsInfections.R')

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# data("cov19mort2021")
# ncast <- nowcasting(data = cov19mort2021, day.effect = FALSE)
# plot(ncast) # Show nowcasted cases
# plot(ncast, type = "delay") # Show contour of delay distribution

# ^ the format of this data is ... ?
# you'have to know the delay distribution first, and
# know the values that haven't been

# see here: https://epilps.com/Nowcasting.html

get_epilps_output <- function(k) {

  # first nowcast the data
  # its a little unclear to see how to do this in
  # https://epilps.com/NowcastingRt.html

  # you can't directly control the # of B-Splines using the nowcasting formulation
  # and its not very clear how to create the data format from
  # aggregated data alone

  # ok so you need to transform this using the nowcast function
  # but its not clear how to do this if you already have aggregated data
  # maybe this isn't a likely scenario but it seems like you'd have to know
  # what the "not yet reported" are to plot this
  # very unclear

  ## reverse the repording PMF to get the amount that will be truncated
  tail_Px <- rev(cumsum(reporting_delay_pmf))
  px_len <- length(tail_Px)

  reports_df$Px = 1
  cases_len <- nrow(reports_df)
  reports_df$Px[(cases_len - px_len + 1):cases_len] <- tail_Px

  ## estimate the expected value
  reports_df$N_EST <- reports_df$N / reports_df$Px

  ## # ok now create the matrix that nowcasting needs
  ## so at the tail, there are limited options
  ## if you are 5 days out, the only way you can be not reported
  ## is if your reporting delay is 5
  ## if you are 4 days out, you can have delays 4 or 5

  ## NOTE: THIS REQUIRES THAT reporting delay has a 0-day
  reports_l <- vector("list", nrow(reports_df))
  set.seed(123)
  for(i in 1:nrow(reports_df)) {
    test_dt <- reports_df$date[i]
    n  <- reports_df$N_EST[i]
    lenPx <- length(reporting_delay_pmf)
    d <- sample(0:(length(reporting_delay_pmf)-1),
                replace = T,
                size = n,
                prob = reporting_delay_pmf)
    reports_l[[i]] <- data.frame(test_dt = test_dt, d = d, t = i)
  }
  linelist <- do.call(rbind, reports_l)
  linelist$report_dt <- linelist$test_dt + linelist$d
  setDT(linelist)

  linelist_agg <- linelist[, .N, by = .(test_dt, d, t, report_dt)][order(t, d)]

  linelist_agg$Reported <- ifelse(linelist_agg$report_dt <=
                                    max(reports_df$date), 'Reported',
                                  'Not yet reported')

  setnames(linelist_agg, c('N', 'test_dt', 'report_dt'),
           c('Cases', 'Date', 'Rep.date'))

  linelist_agg <- data.frame(linelist_agg[, c('t', 'd', 'Date', 'Rep.date',
                                   'Cases', 'Reported')][order(t, d)])
  head(linelist_agg)

  # using the MCMC version
  LPSfit <- nowcastingR(data = linelist_agg,
                        method = k,
                   si = serial_interval_pmf)

  # create dataframe
  # shift backwards by the seeding_time, aka the sum
  # of the means of the delay distributions
  plot_epilps <- data.frame(
    model = paste0("k = 40"),
    ###
    date = LPSfit$Rnow$Time, # - seeding_time,
    ###
    Rt = LPSfit$Rnow$Rq0.50,
    Rt_lb = LPSfit$Rnow$Rq0.025,
    Rt_ub = LPSfit$Rnow$Rq0.975
  )
  return(plot_epilps)
}

df1 <- get_epilps_output("M3")
epilps_df <- rbind(df1)

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

ggsave('img/BSplines.png', width = 6.5, height = 1.5)
saveRDS(plot_rt1, "img/BSplines.RDS")

