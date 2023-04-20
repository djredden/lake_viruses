rm(list = ls())
library(tidyverse)
library(ggtext)

data_cleaned <- read_csv("cleaned_data/data_cleaned.csv")

time_series <- data_cleaned %>% 
  #filter(target != "Measles") %>% 
  mutate(sample_date = lubridate::ymd(sample_date),
         gu_sample = if_else(is.na(gu_sample), 1e-4, gu_sample)) %>% 
  ggplot(aes(sample_date, log(gu_sample), color = sample_id, shape = adsorbent)) +
  geom_point(size = 3, aes(alpha = detect)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  facet_wrap(vars(target)) +
  ggsci::scale_color_jco() +
  theme_bw() +
  labs(y = "Log Gene Target Concentration <br> (GU g<sup>-1</sup>)",
       x = NULL,
       shape = "Sample Type") +
  theme(axis.title.y = element_markdown()) +
  guides(alpha = "none")

shape_legend <- cowplot::get_legend(time_series + guides(color = "none") +
                                      theme(legend.background = element_rect(fill = "white", color = "white")))

plot_combined <- cowplot::plot_grid(time_series +
                                      guides(shape = "none") +
                                      theme(legend.position = "bottom") +
                                      labs(color = NULL),
                                    shape_legend,
                                    ncol = 2, 
                                    rel_widths = c(.8, .2)) 

cowplot::ggdraw(plot_combined) + 
  theme(plot.background = element_rect(fill="white", color = NA))


ggsave("output/time_series.png", width = 8, height = 7)



