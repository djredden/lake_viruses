rm(list = ls())
library(tidyverse)
library(scales)
library(ggtext)

theme_set(theme_bw() +
            theme(
              axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 22),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 20),
              strip.text = element_text(size = 22, face = "bold"),
              panel.spacing.x = unit(2, "lines"),
              panel.spacing.y = unit(1, "lines")
            ))

data_cleaned <- read_csv("cleaned_data/data_cleaned.csv")

# --------------------------------------------------- Detection frequencies ---------------------------------------------------

det_freq <- data_cleaned %>% 
  filter(adsorbent == "GAC" & extraction == "tween_bead" |
                    adsorbent == "Grab") %>% 
  group_by(target, adsorbent) %>%
  summarise(
    n = n(),
    n_det = sum(detect),
    freq_detect = 100*mean(detect),
    freq_detect = round(freq_detect, 1)) %>% 
  ungroup() 

write_csv(det_freq, "output/detection_frequencies.csv")

# -------------------------------------------------------- Time series plot --------------------------------------------------------

time_series <- data_cleaned %>% 
  filter(detect,
         adsorbent == "GAC",
         extraction == "tween_bead") %>% 
  left_join(., det_freq, by = c("target","adsorbent")) %>% 
  ggplot(aes(sample_date, gu_total, color = sample_id)) +
  geom_point(size = 3) +
  facet_wrap(vars(target)) +
  scale_color_manual(values = c("#0073C2FF", "darkgrey")) +
  
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +

  #scale_y_continuous(trans ='log10', labels = math_format()) +

  labs(y = "Total Gene Copies (GC)",
       x = NULL,
       color = NULL) +
  theme(axis.title.y = element_markdown(margin = margin(t = 0, r = 12, b = 0, l = 0)),
        strip.text = element_markdown(lineheight = 1.2,
                                      halign = 0.5,
                                      padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
        legend.position = "bottom",
        plot.margin = margin(10, 20, 10, 10)) 


ggsave("output/time_series.png", height = 8, width = 10, dpi = 300)

#  --------------------------------------------------- Barplot of detections ---------------------------------------------------

data_cleaned %>%  
  filter(adsorbent == "GAC" & extraction == "tween_bead" |
           adsorbent == "Grab") %>% 
  group_by(target, adsorbent) %>% 
  summarise(freq_detect = 100*mean(detect)) %>% 
  ungroup() %>% 
  mutate(target = fct_reorder(target, freq_detect, max, .desc = TRUE)) %>% 
  # Change zeros to small negative value for plotting purposes
  mutate(freq_detect = if_else(freq_detect == 0, -0.5, freq_detect)) %>% 
  ggplot(aes(target, freq_detect, fill = adsorbent)) +
  geom_col(position = position_dodge(),
           width = .8) +
  ggsci::scale_fill_jco() +
  theme(legend.position = "bottom",
        axis.text.x = element_markdown()) +
  labs(x = NULL,
       y = expression(paste("Positive Gene Target Detections (%)")),
       fill = NULL) +
  coord_cartesian(ylim = c(0,100)) +
  theme(legend.position = c(0.87, 0.9))

ggsave(here::here("output/detection_frequency.png"), height = 8, width = 12, dpi = 300)


