rm(list = ls())
library(tidyverse)

# Read qPCR data and put in long format
raw_data <- read_csv("raw_data/data_raw_full.csv") %>% 
  mutate(across(
    .cols = contains("_cq"),
    .fns = as.numeric
  )) %>% 
  pivot_longer(cols = contains("_cq"),
               names_to = "target", 
               values_to = "cq") %>% 
  mutate(target = str_remove(target, "_cq"),
         cq = if_else(cq == 0, NA_real_, cq),
         detect = if_else(is.na(cq) | cq >= 38, FALSE, TRUE)) 


# Read in standard curve results
std_curve_data <- read_csv("raw_data/std_curves.csv")


# Different extraction methods used
extraction_types <- raw_data %>% 
  distinct(extraction) %>% 
  filter(!is.na(extraction)) %>% 
  pull()

  
# Combine raw data and curve results
data_cleaned <- left_join(raw_data, std_curve_data, by = "target") %>% 
  # apply standard curve to calculate gene concentration in lysate (3uL template used per rxn)
  mutate(gu_lysate = 10^((cq - y_int) / slope)/3) %>% 
  select(-c(r2:efficiency)) %>% 
  # combine with processing metadata
  mutate(gu_total = case_when(
    adsorbent == "gac" & is.na(extraction) ~ gu_lysate * lysis_vol_ul,
    adsorbent == "grab" & is.na(extraction) ~ gu_lysate * lysis_vol_ul,
    adsorbent == "gac" & extraction %in% extraction_types ~ (gu_lysate * (eluate_vol_ul + lysis_vol_ul) / eluate_vol_ul) * eluent_vol_ml * 1000,
    adsorbent == "grab" & extraction %in% extraction_types ~ gu_lysate * eluent_vol_ml * 1000,
    TRUE ~ NA
  ),
  gu_sample = case_when(
    adsorbent == "gac" & is.na(extraction) ~ gu_lysate * lysis_vol_ul / ads_mass_mg,
    adsorbent == "grab" & is.na(extraction) ~ gu_lysate * lysis_vol_ul / volume_ml,
    adsorbent == "gac" & extraction %in% extraction_types ~ (gu_lysate * (eluate_vol_ul + lysis_vol_ul) / eluate_vol_ul) * eluent_vol_ml * 1000  / ads_mass_mg,
    adsorbent == "grab" & extraction %in% extraction_types ~ gu_lysate * eluent_vol_ml * 1000 / ads_mass_mg,
    TRUE ~ NA
  )
  ) %>% 
  mutate(target = factor(target, levels = c("mev", "infa", "rsv", "sars", "adv", "nv", "env", "rv")),
         target = case_when(
           target == "adv" ~ "Adenovirus",
           target == "rv" ~ "Rotavirus",
           target == "sars" ~ "SARS-CoV-2",
           target == "infa" ~ "Influenza A",
           target == "nv" ~ "Norovirus",
           target == "env" ~ "Enterovirus",
           target == "rsv" ~ "RSV",
           target == "mev" ~ "Measles"
         )) %>% 
  mutate(sample_id = factor(sample_id, levels = c("banook1", "banook2")),
         sample_id = if_else(sample_id == "banook1", "Site 1", "Site 2")) %>% 
  mutate(adsorbent = if_else(adsorbent == "gac", "GAC", "Grab")) %>% 
  mutate(sample_date = lubridate::ymd(sample_date)) 
  
write_csv(data_cleaned,  "cleaned_data/data_cleaned.csv")

