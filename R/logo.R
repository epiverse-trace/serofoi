
library(tidyverse)
rm(list=ls())

dat <- readRDS("data/dat_plot.RDS")


ggplot(dat) +
  geom_ribbon(aes( x= age, ymin = predicted_prev_lower, ymax = predicted_prev_upper),
              fill = '#41b6c4', alpha = 0.8) +
  geom_line(aes( x= age, y = predicted_prev), colour = '#0570b0', size = 5) +
  geom_errorbar(aes(age, ymin = observed_prev_lower, ymax = observed_prev_upper), width = 0,
                colour = "#41ab5d", size = 3) +
  geom_point(aes(age, observed_prev), size = 13, fill = '#f768a1', colour = 'black', shape = 21) +
  theme_void() +
  coord_cartesian(xlim = c(0, 65), ylim = c(0,1.05)) +
  theme(legend.position = 'none') +
  ylab ('') + xlab("")

dev.off()


