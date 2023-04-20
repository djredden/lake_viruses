rm(list = ls())
library(tidyverse)
library(ggtext)

data_cleaned <- read_csv("cleaned_data/data_cleaned.csv")

data_cleaned %>% 
  group_by(target, sample_id, adsorbent) %>% 
  summarise(freq_detect = 100*mean(detect)) %>% 
  ungroup() %>% 
  mutate(freq_detect = if_else(freq_detect == 0, -0.5, freq_detect)) %>% 
  
  ggplot(aes(sample_id, freq_detect, fill = adsorbent, group = target)) +
  geom_col(position = position_dodge2(preserve = "single"),
           width = 1) +
  ggsci::scale_fill_jco() +
  facet_wrap(vars(target)) +
  theme(legend.position = "bottom",
        axis.text.x = element_markdown()) +
  labs(x = NULL,
       y = expression(paste("Positive Gene Target Detections (%)")),
       fill = NULL) +
  theme_bw() 

ggsave(here::here("output/detection_frequency.png"), height = 8, width = 8, dpi = 300)

