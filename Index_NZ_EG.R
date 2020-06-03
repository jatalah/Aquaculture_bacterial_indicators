library(tidyverse)
library(readxl)

dat_med_log <-  
  read_csv('data/bacteria.csv',col_types = cols()) %>%
  select(-c(2:12)) %>% 
  gather(ASV, abund, -SampleID)


seq_codes <- read_excel('data/Barcode_ASVnumber.xlsx') %>% 
  rename(Sequence = ASV_Barcode) %>% 
  select(ASV,Sequence)
  
nz_seq_codes <- read_csv('data/EG_groups_bact_NZ_taxo_wSeq.csv') %>% 
  select(Sequence,EG)

EG_NZ_Med <- 
  inner_join(seq_codes,nz_seq_codes, by = "Sequence") %>% 
  select(-Sequence)


# EG proportion calculations --------------
eg_prop_nz_eg <- 
  dat_med_log %>%
  inner_join(EG_NZ_Med, by = 'ASV') %>%
  group_by(SampleID) %>%
  mutate(N = sum(abund)) %>% 
  group_by(SampleID, EG) %>%
  mutate(n = sum(abund) / N) %>%
  summarise(n = first(n) * 100) %>%
  spread(EG, n)


# Combine all indices with env----------------
env <- read_csv('data/env_bact.csv')

AMBI_weights <- lm(TFS ~  0 + I + II + III + IV + V, data = EG_prop_all_nz)

EG_prop_all_nz <- left_join(eg_prop_nz_eg, env, by = 'SampleID') %>% 
  mutate(
    bMBI = (
      0.017 * I +
        0.030  * II + 
        0.025 * III + 
        0.039 * IV + 
        0.066 * V
    ),
    bMBI_local = (
      AMBI_weights$coefficients[1] * I +
        AMBI_weights$coefficients[2] * II + 
        AMBI_weights$coefficients[3] * III + 
        AMBI_weights$coefficients[4] * IV + 
        AMBI_weights$coefficients[5] * V
    )/100
  ) %>% 
  ungroup() %>% 
  mutate(Distance = fct_relevel(Distance, "0m", "50m", "100m", "250m", "500m", "1000m","2000m", "Control" ))

# Plot euk-AMBI vs macrofaunal AMBI-----------
ggplot(EG_prop_all_nz, aes(Distance, bMBI_local, color = Farm)) +
  geom_boxplot() +
  labs(y = 'b-MBI', x =  'Distance')

ggplot(EG_prop_all_nz, aes(TFS, bMBI)) +
  geom_point(size = 4, aes(color  = Distance, shape = Farm)) +
  stat_smooth(method = "lm", alpha = .2) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = ..rr.label..),
    parse = TRUE,
    size = 4
  ) +
  scale_color_brewer(palette = 'Spectral') +
  labs(y = 'b-MBI', x =  'Total free sulphide (TFS)')

ggplot(EG_prop_all_nz, aes(TFS, bMBI_local)) +
  geom_point(size = 4, aes(color  = Distance, shape = Farm)) +
  stat_smooth(method = "lm", alpha = .2) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = ..rr.label..),
    parse = TRUE,
    size = 4
  ) +
  scale_color_brewer(palette = 'Spectral') +
  labs(y = 'b-MBI', x =  'Total free sulphide (TFS)')

