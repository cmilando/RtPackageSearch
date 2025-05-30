library(patchwork)

pA <- readRDS("img/Infections.RDS") +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pB <- readRDS("img/FixedSlidingWindow.RDS") +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pC <- readRDS("img/RandomWalk.RDS") +
  theme(
        panel.grid.major = element_line(color = grey(0.9),
                                        linewidth = 0.1))

pA / pB / pC

##
ggsave('img/Fig3.png', width = 6.5, height = 1.5 * 3)

##
ggsave('img/Fig3.png', width = 8, height = 5)

## probabily somewhere between these two is what i want
