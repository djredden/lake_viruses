rm(list = ls())
library(tidyverse)
library(janitor)
library(ggtext)
library(kableExtra)


water_chem <- read_csv("raw_data/water_chemistry.csv") %>% 
  clean_names() %>% 
  select(-sample_type, -ysi_id, -feem_file_id, -starts_with(c("tn", "x")),
         -nitrate, -nitrite, -phosphate, -sulfate, -chloride
  ) %>% 
  mutate(across(
    .cols = -c(sample_id, date),
    .fns = as.numeric
  ))

data_long <- water_chem %>% 
  # Take means of replicates
  rowwise() %>% 
  mutate(
    turbidity = mean(c(turbidity_1, turbidity_2), na.rm = TRUE),
    colour = mean(c(colour_1, colour_2), na.rm = TRUE),
    uv_254 = mean(c(uv_254_1, uv_254_2), na.rm = TRUE)
  ) %>% 
  # Remove original replicates
  select(-matches("_(1|2)$")) %>% 
  # Put in long format
  pivot_longer(cols = -c(sample_id, date),
               names_to = "param",
               values_to = "value") %>% 
  mutate(param = case_when(
    param == "al" ~ "Total Aluminum (&mu;g L<sup>-1</sup>)",
    param == "colour" ~ "Colour (Pt-Co)",
    param == "doc" ~ "DOC (mg L<sup>-1</sup>)",
    param == "fe" ~ "Total Iron (&mu;g L<sup>-1</sup>)",
    param == "field_cond" ~ "Conductivity (&mu;s cm<sup>-1</sup>)",
    param == "field_do" ~ "DO (mg L<sup>-1</sup>)",
    param == "field_ph" ~ "pH",
    param == "field_temp" ~ "Temperature (&deg;C)",
    param == "tds"  ~ "TDS (mg L<sup>-1</sup>)",
    param == "toc" ~ "TOC (mg L<sup>-1</sup>)",
    param == "p" ~ "Total P (&mu;g L<sup>-1</sup>)",
    param == "turbidity" ~ "Turbidity (NTU)",
    param == "uv_254" ~ "UV<sub>254</sub> (cm<sup>-1</sup>)"),
  sample_id = fct_recode(sample_id,
                         "Site 1" = "banook1",
                         "Site 2" = "banook2")
  ) 
  
  

# --------------------------------------------- Plot of all water chemistry ---------------------------------------------

data_long %>% 
  ggplot(aes(x = sample_id, y = value)) +
  geom_boxplot() +
  facet_wrap(vars(param), scales = "free") +
  theme_bw() +
  labs(y = NULL,
       x = NULL) +
  theme(strip.text = element_markdown())

ggsave("output/water_chem_by_site.png", dpi = 300,
       height = 8, width = 10)


# ----------------------------------------------------------- Summary Table -----------------------------------------------------------

summary_long <- water_chem %>% 
  # Take means of replicates
  rowwise() %>% 
  mutate(
    turbidity = mean(c(turbidity_1, turbidity_2), na.rm = TRUE),
    colour = mean(c(colour_1, colour_2), na.rm = TRUE),
    uv_254 = mean(c(uv_254_1, uv_254_2), na.rm = TRUE)
  ) %>% 
  # Remove original replicates
  select(-matches("_(1|2)$")) %>% 
  # Put in long format
  pivot_longer(cols = -c(sample_id, date),
               names_to = "param",
               values_to = "value") %>% 
  mutate(param = case_when(
    param == "al" ~ "Total Aluminum (µg/L)",
    param == "colour" ~ "Colour (Pt-Co)",
    param == "doc" ~ "DOC (mg/L)",
    param == "fe" ~ "Total Iron (µg/L)",
    param == "field_cond" ~ "Conductivity (µs/cm)",
    param == "field_do" ~ "DO (mg/L)",
    param == "field_ph" ~ "pH",
    param == "field_temp" ~ "Temperature (C)",
    param == "tds"  ~ "TDS (mg/L)",
    param == "toc" ~ "TOC (mg/L)",
    param == "p" ~ "Total P (µg/L)",
    param == "turbidity" ~ "Turbidity (NTU)",
    param == "uv_254" ~ "UV254 (/cm)"),
    sample_id = fct_recode(sample_id,
                           "Site 1" = "banook1",
                           "Site 2" = "banook2")
  ) %>% 
  # Summarise
  group_by(param) %>% 
  summarise(
    Min = min(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE)) %>% 
  mutate(
    across(where(is.numeric),
           ~round(.x, 2)) 
  ) %>% 
  rename("Parameter" = param)

write_csv(summary_long, "output/water_chemistry_table.csv")

