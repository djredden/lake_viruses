rm(list = ls())
library(tidyverse)
library(ggtext)

theme_set(theme_bw() +
            theme(
              axis.text = element_text(size = 20, color = "black"),
              axis.title = element_text(size = 22),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 20),
              strip.text = element_text(size = 22, face = "bold"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1, "lines")
            ))

data_cleaned <- read_csv("cleaned_data/data_cleaned.csv")  %>% 
  # Other methods were used to extract NAs from GAC. Only one is included in this work
  filter(adsorbent == "GAC" & extraction == "tween_bead" |
           adsorbent == "Grab") 

# ----------------------------------------------- Summary stats of sampling -----------------------------------------------

n_events <- data_cleaned %>% 
  group_by(sample_id, adsorbent, target) %>% 
  summarise(
    n= n()
  ) %>% 
  distinct(sample_id, adsorbent, n)

write_csv(n_events, "output/number_sample_events.csv")

# ------------------------------------------------------- Wilcoxon matched pairs -------------------------------------------------------

# Using wilcox.test function
paired_wilcoxon <- data_cleaned %>%
  filter(adsorbent == "GAC") %>% 
  select(target, sample_date, sample_id, gu_lysate) %>% 
  pivot_wider(id_cols = c(target, sample_date),
              names_from = sample_id,
              values_from = gu_lysate) %>% 
  mutate(
    across(
      where(is.numeric),
      ~replace_na(., 0))
  ) %>% 
  group_by(target) %>% 
  nest() %>% 
  mutate(model = map(data, ~wilcox.test(.$`Banook 1`, .$`Banook 2`, paired = TRUE, conf.int = TRUE)),
         tidy = map(model, broom::glance)) %>% 
  unnest(tidy) %>% 
  select(-c(data, model))

write_csv(paired_wilcoxon, "output/difference_estimates.csv")

# ------------------------------------------------- Boxplots of differences -------------------------------------------------

difference_plot <- data_cleaned %>%
  filter(adsorbent == "GAC") %>% 
  select(target, sample_date, sample_id, gu_lysate) %>% 
  ggplot(aes(sample_id, gu_lysate)) +
  geom_boxplot() +
  facet_wrap(vars(target), scales = "free_y")
