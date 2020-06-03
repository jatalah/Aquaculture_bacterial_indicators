library(tidyverse)
indices_all <- read_csv('data/indices_all.csv')

# fit lms linear, sqrt and log-----------
indices_lms <-
  indices_all %>%
  gather(index, value, bMBI_JA:bMBI_borja) %>%
  group_by(index) %>%
  nest() %>%
  mutate(
    lm = map(data, ~ lm(value ~ TFS, data = .x)),
    lm_sqrt = map(data, ~ lm(value ~ sqrt(TFS), data = .x)),
    lm_log = map(data, ~ lm(value ~ log(TFS), data = .x)),
    lm_table = map(lm, tidy),
    lm_glance = map(lm, glance),
    lm_sqrt_table = map(lm_sqrt, tidy),
    lm_sqrt_glance = map(lm_sqrt, glance),
    lm_log_table = map(lm_log, tidy),
    lm_log_glance = map(lm_log, glance)
  )

# save summary tables---------
indices_lms %>%
  select(lm_table:lm_log_glance) %>%
  unnest() %>%
  ungroup() %>%
  write_csv('outputs/index_lm_tables.csv')

# save the index data with Keely formula----
indices_all %>%
  select(SampleID, bMBI_keely) %>%
  rename(AMBI_bact = "bMBI_keely") %>%
  write_csv('outputs/bMBI_data.csv')