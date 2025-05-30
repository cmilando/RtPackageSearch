# Create incidence time-series for reports and infections

#
library(ggplot2)
library(ggpubr)
library(data.table)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# SIMULATE INFECTIONS
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# create some infections data
create_linelist_data <- function() {
  set.seed(123)
  ## peak 1 - longer
  days      <- 1:100
  amplitude <- 0.55
  period    <- 100
  frequency <- 4
  sine_wave <- amplitude * sin(2 * pi * frequency * days / period)
  ## peak 2 - sharp
  days      <- 1:50
  amplitude <- 0.65
  period    <- 50
  frequency <- 4
  sine_wave2 <- amplitude * sin(2 * pi * frequency * days / period)
  ## peak 3 - longer
  new_ts    <- rev(c(sine_wave[1:25],
                     rep(sine_wave2[23], 10),
                     sine_wave2[12:23],
                     rep(sine_wave2[23], 10),
                     sine_wave[46:70],
                     rep(sine_wave[70], 10)))
  lambda    <- exp(new_ts[1:82]) * 1000 - 300

  cases     <- sapply(lambda, function(l) rpois(1, l + exp(rnorm(1, 0, 3))))
  ## a little splice
  cases[5:6] <- cases[1:2]
  cases[68:71] <- c(300, 312, 305, 330)
  cases[76] <- 270
  ##
  dt        <- seq.Date(from = '2020-03-01', length.out = length(cases))
  df        <- data.frame(date = dt, N = cases)
  df_l      <- lapply(1:nrow(df), function(i) rep(df$date[i], df$N[i]))
  df_c      <- do.call(c, df_l)
  df_r      <- data.table(data.frame(date = df_c))
  return(df_r)
}

linelist <- create_linelist_data()

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
# this combines infection -> onset -> test
# has to start with 0 !
infect_to_test_pmf <- c(
  0, 0.0001, 0.0030, 0.1422, 0.2714, 0.2664,
  0.1832, 0.0846, 0.0340, 0.0136, 0.001, 0.0005
)

# 2) time from taking a test to it getting into a state database
reporting_delay_pmf <- c(
  0.3786, 0.3724, 0.1662, 0.0622, 0.0166, 0.004
)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# CONVOLVE TO INFECTIONS FROM REPORTS
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

set.seed(123)

total_delay <-
  # take a draw from this distribution
  sample(0:(length(infect_to_test_pmf) - 1),
         size = nrow(linelist),
         prob = infect_to_test_pmf, replace = TRUE) +
  # and this one
  sample(0:(length(reporting_delay_pmf) - 1),
         size = nrow(linelist),
         prob = reporting_delay_pmf, replace = TRUE) +
  # some random noise
  sample(c(-2:2), size = nrow(linelist), replace = T)


get_mean <- function(pmf, start) {
  stopifnot(start %in% c(0, 1))
  weighted.mean(x = start:(length(pmf) - 1 + start),
                w = pmf)
}

## some error with seeding time here but
## i'll calculate it outside
w1 <- get_mean(generation_interval_pmf, 0)
w2 <- get_mean(infect_to_test_pmf, 0)
w3 <- get_mean(reporting_delay_pmf, 1)
seeding_time <- round(w1 + w2 + w3)

linelist$report_date <- linelist$date + total_delay

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# AGGREGATE TO DAILY
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

#
infections_df <- linelist[, .N, by = date][order(date)]
reports_df <- linelist[, .N, by = report_date][order(report_date)]

#
setnames(reports_df, 'report_date', 'date')

#
reports_df$type <- 'Observed\nreports'
infections_df$type <- 'Unobserved\ninfections'

#
burn_in <- 6
reports_df <- subset(reports_df,
                     date >= min(infections_df$date) + burn_in &
                     date <= max(infections_df$date))

infections_df <- subset(infections_df,
                     date >= min(infections_df$date) + burn_in)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# ADD RIGHT TRUNCATION TO DATA BASED ON REPORTING DELAY
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

tail_Px <- rev(cumsum(reporting_delay_pmf))
px_len <- length(tail_Px)

reports_df$Px = 1
cases_len <- nrow(reports_df)

reports_df$Px[(cases_len - px_len + 1):cases_len] <- tail_Px

reports_df$N      <- reports_df$N * reports_df$Px
reports_df$N      <- as.integer(reports_df$N)
reports_df$Px <- NULL

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////


##
report_ymax = 2000
first_day <- min(reports_df$date)
last_day <- max(reports_df$date)
nowcast_start    = last_day - seeding_time
forecast_window  = last_day + 15

##

plot_inf <- ggplot(rbind(infections_df, reports_df)) +
  theme_classic2() +
  # ##
  geom_line(aes(x = date, y = N, color = type,
                linetype = type),
            linewidth = 0.5, show.legend = T) +
  scale_linetype_manual(name = NULL, values = c('solid', '11')) +

  geom_point(aes(x = date, y = N, color = type, shape = type),
             size = 1, fill = 'white',
             show.legend = T) +
  scale_shape_manual(name = NULL, values = c(21, NA)) +
  scale_color_manual(name = NULL, values = c(grey(0.75), 'pink')) +
  ##
  coord_cartesian(xlim = c(first_day,
                           forecast_window),
                  ylim = c(0, report_ymax+100), expand = F) +
  ylab("Daily cases") +
  xlab(NULL) +
  ##
  scale_x_date(breaks = '2 week',
               date_minor_breaks = "1 weeks",
               date_labels = "%b %d") +
  ##
  annotate('text', x = first_day + 3,
           color = 'blue',
           y = report_ymax, label = 'A.',
           fontface= 'bold',
           size = 5) +
  annotate('text', x = first_day + 35,
           y = report_ymax, label = 'Historical period',
           size = 2) +
  annotate('text', x = nowcast_start + 6,
           y = report_ymax, label = 'Nowcasting', size = 2) +
  annotate('text', x = last_day + 6,
         y = report_ymax, label = 'Forecasting', size = 2) +
  annotate('text', x = first_day + 11,
           color = scales::muted("red"),
           y = 1400, label = 'Peak 1', size = 2.5) +
  annotate('text', x = first_day + 41,
           color = scales::muted("red"),
           y = 1400, label = 'Peak 2', size = 2.5) +
  annotate('segment',
           arrow = arrow(type = "closed", length = unit(0.03, "npc")),
           color = grey(0.5),
           linewidth = 0.5,
           x = nowcast_start + 7.5,
           xend = nowcast_start + 9.5,
           y = 140,
           yend = 220) +
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
           y = 1500,
           angle = 90,
           label = 'PRESENT') +
  theme(legend.key.spacing.y = unit(.25, 'cm'))

plot_inf
ggsave('img/Infections.png', width = 6.5, height = 2)
saveRDS(plot_inf, "img/Infections.RDS")
