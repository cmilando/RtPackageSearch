library(patchwork)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# RUN ALL
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////
rm(list = ls()); gc();
source('01_ReportsInfections.R'); rm(list = ls()); gc();
source('02_FixedSlidingWindows.R'); rm(list = ls()); gc();
source('03_RandomWalk.R'); rm(list = ls()); gc();
source('04_Filtering.R'); rm(list = ls()); gc();
source('05_BSplines.R'); rm(list = ls()); gc();
source('06_GaussianProcess.R'); rm(list = ls()); gc();

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

source('02_FixedSlidingWindows.R')

pA <- readRDS("img/Infections.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           hjust = 0,
           y = report_ymax,
           label = 'A. Observed reports',
           fontface= 'bold',
           size = 3) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.spacing.y = unit(.1, 'cm'),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pB <- readRDS("img/GaussianProcess.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           y = rt_max,
           hjust = 0,
           label = 'B. Gaussian Process (EpiNow2, epinowcast)',
           fontface= 'bold',
           size = 3) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pC <- readRDS("img/RandomWalk.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           y = rt_max,
           hjust = 0,
           label = 'C. Random walk (EpiNow2, epinowcast)',
           fontface= 'bold',
           size = 3) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pD <- readRDS("img/FixedSlidingWindow.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           y = rt_max,
           hjust = 0,
           label = 'D. Sliding window (EpiEstim)',
           fontface= 'bold',
           size = 3) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pE <- readRDS("img/Filtering.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           y = rt_max,
           hjust = 0,
           label = 'E. Trend-filtering (RtEstim)',
           fontface= 'bold',
           size = 3) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = grey(0.9),
                                    linewidth = 0.1))

pF <- readRDS("img/BSplines.RDS") +
  annotate('text', x = first_day + 1,
           color = 'blue',
           y = rt_max,
           hjust = 0,
           label = 'F. P-splines (EpiLPS)',
           fontface= 'bold',
           size = 3) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = grey(0.9),
                                    linewidth = 0.1))

pA / pB / pC / pD / pE / pF

##
ggsave('img/Fig3.png', width = 8, height = 8)

