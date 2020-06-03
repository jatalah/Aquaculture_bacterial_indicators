# Calculate the bacterial AMBI------------------
library(tidyverse)

# read AVS read abundance and eco-groups data-----------------------------------
data <-
  read_csv('data/bacteria.csv', col_types = cols()) %>%
  select(-c(2:12))

EG_groups <- read_csv('outputs/EG_groups_MED.csv')

# EG proportion calculations --------------
eg_prop <-
  data %>%
  gather(ASV, abund,-SampleID) %>%
  inner_join(EG_groups, by = 'ASV') %>%
  group_by(SampleID) %>%
  mutate(N = sum(abund)) %>%
  group_by(SampleID, EG) %>%
  mutate(n = sum(abund) / N) %>%
  summarise(n = first(n) * 100) %>%
  spread(EG, n)


# Combine all indices with env----------------
env <- read_csv('data/env_bact.csv')

EG_prop_all <- left_join(eg_prop, env, by = 'SampleID')

# Calculate optimal weights for AMBI using linear regression on the train dataset-------
AMBI_weights <-
  lm(TFS ~   0 + I + II + III + IV + V, data = EG_prop_all)
AMBI_weights
# lm(TFS ~  I + III +  V, data = EG_prop_all)

indices_all <-
  EG_prop_all %>%
  mutate(
    bMBI_JA = (
      AMBI_weights$coefficients[1] * I +
        AMBI_weights$coefficients[2] * II +
        AMBI_weights$coefficients[3] * III +
        AMBI_weights$coefficients[4] * IV +
        AMBI_weights$coefficients[5] * V
    ) / 100,
    bMBI_keely = (1.5 * II + 3 * III + 5 * IV + 12 * V) / 100,
    
    bMBI_borja = (1.5 * II + 3 * III + 4.5 * IV + 6 * V) / 100
    
  ) %>%
  ungroup() %>%
  mutate(Distance = str_replace(Distance, "m"," m"),
         Distance = fct_relevel(
           Distance,
           "0 m",
           "50 m",
           "100 m",
           "250 m",
           "500 m",
           "1000 m",
           "2000 m",
           "Control"
         )) %>% 
  mutate(Farm = fct_recode(Farm, `Farm 1` = "Ext", `Farm 2` = "Int")) 

write_csv(indices_all, 'data/indices_all.csv')
