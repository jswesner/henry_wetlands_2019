library(brms)
library(tidyverse)
library(tidybayes)

brm_bug_selenium = readRDS("models/brm_bug_selenium.rds")
brm_water = readRDS("models/brm_water.rds")
sites_lat_lon = readRDS(file = "data/sites_lat_lon.rds")

bug_posteriors = brm_bug_selenium$data %>% 
  distinct(taxon, trt, mon_yr, site) %>%
  mutate(month = case_when(grepl("July", mon_yr) ~ "July",
                           grepl("June", mon_yr) ~ "June",
                           TRUE ~ "May"),
         year = parse_number(mon_yr)) %>% 
  mutate(stage = "adults",
         'site:rep' = "new") %>% 
  add_epred_draws(brm_bug_selenium, re_formula = NULL, allow_new_levels = T) %>% 
  group_by(.draw, trt, site, year, month) %>% 
  reframe(.epred = mean(.epred)) %>% 
  mutate(measure = "se_ugg")

water_posteriors = brm_water$data %>% 
  distinct(date2, site, trt) %>%
  mutate(month = case_when(grepl("July", date2) ~ "July",
                           grepl("June", date2) ~ "June",
                           TRUE ~ "May"),
         year = parse_number(date2)) %>% 
  add_epred_draws(brm_water, re_formula = NA) %>% 
  mutate(measure = "water_ugl")
  

bug_water_posteriors = bind_rows(bug_posteriors, water_posteriors) %>% 
  mutate(trt = str_to_lower(trt),
         site = str_to_lower(site))


aquasync_summary = bug_water_posteriors %>% 
  group_by(site, trt, month, year, measure) %>% 
  reframe(mean = mean(.epred),
          sd = sd(.epred)) %>% 
  mutate(mean_sd = paste0(mean, "_", sd)) %>% 
  select(-mean, -sd) %>% 
  pivot_wider(names_from = measure, values_from = mean_sd) %>% 
  filter(!is.na(se_ugg)) %>% 
  filter(!is.na(water_ugl)) %>% 
  separate(water_ugl, into = c("water_ugl_mean", "water_ugl_sd"), sep = "_") %>% 
  separate(se_ugg, into = c("bug_ugg_mean", "bug_ugg_sd"), sep = "_") %>% 
  mutate(water_ugl_mean = as.numeric(water_ugl_mean),
         water_ugl_sd = as.numeric(water_ugl_sd),
         bug_ugg_mean = as.numeric(bug_ugg_mean),
         bug_ugg_sd = as.numeric(bug_ugg_sd)) %>% 
  mutate(sample_name = paste(site, trt, month, year, sep = "_")) %>% 
  select(sample_name, water_ugl_mean, water_ugl_sd, 
         bug_ugg_mean, bug_ugg_sd, everything()) %>% 
  left_join(sites_lat_lon %>% mutate(site = str_to_lower(site)) %>% select(site, lat, long))


write_csv(aquasync_summary, file = "data/aquasync_summary.csv")


