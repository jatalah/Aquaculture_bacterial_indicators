# load libraries------------
library(tidyverse)
library(splines)
library(quantreg)
library(knitr)
library(pracma)# for finding peaks
library(ggpubr)
library(vegan)
library(phyloseq)
library(ggpmisc)

source('theme_javier.R')
theme_set(theme_javier())

#read data --------------
data_all <- read_csv('data/bacteria.csv',col_types = cols())

metadata <- data_all %>% select(1:12)

bact_rare <- 
  data_all %>% 
  column_to_rownames('SampleID') %>% 
  select(-c(1:11)) %>% 
  otu_table(taxa_are_rows = F) %>% 
  rarefy_even_depth(sample.size = 10000, rngseed = T)

ntaxa(bact_rare)

bact_rare_df <-
  as(otu_table(bact_rare), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata, ., by = 'SampleID')


head(names(bact_rare_df), 20)


# select 2000 most abundant ASV------------
abund_ASV <-  
  bact_rare_df %>% 
  gather(ASV, abund, ASV_001:ncol(.)) %>% 
  group_by(ASV) %>% 
  summarise_at(vars(abund), list(sum), na.rm =T) %>%
  top_n(2000) %>% 
  dplyr::select(ASV) %>% 
  as.list()

bact_top <-
  bact_rare_df %>%
  select(Station, Distance, Farm, TFS, abund_ASV$ASV)

# convert into long format
bact_top_long <- 
  bact_top %>% 
  gather(ASV,abund, abund_ASV$ASV)

# create plot with quantile regression and splines with 3 df--------------
p2 <- 
  ggplot(bact_top_long, aes(x = TFS, y = abund, color = Farm)) +
  geom_point(size = 3, alpha = .3) +
  labs(x = 'Total free sulphides', y = 'Number of reads') +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ASV, scales = 'free')


# get splines data from plot------------------------
gb <- ggplot_build(p2)

# obtaine ES values where abundance peak by Farm-------
gb1 <- 
  gb$data[[2]] %>% 
  data.frame() %>% 
  group_by(group,PANEL) %>% 
  nest() %>% 
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(data, ~ findpeaks(
      .$y,
      nups = 0,
      ndowns = 0,
      minpeakheight = max(.$y) / 2.4
    )),
    n_peaks = map_dbl(find_peaks, nrow)
  )


plot_id <- 
  gb$plot[[1]] %>% 
  distinct(Farm,ASV) %>% 
  arrange(ASV)


# get the number of peaks per AVS and farm ------
n_peaks <- 
  gb1 %>% 
  select(n_peaks) %>% 
  unnest() %>% 
  bind_cols(plot_id)


# merge difference peaks by Farm  and the number of peaks-----
peak_dat <-
  gb1 %>%
  select(peak) %>%
  unnest() %>%
  bind_cols(plot_id) %>%
  spread(Farm, peak) %>%
  mutate(diff = abs(Ext - Int)) %>%
  gather(Farm, Peak, Ext:Int) %>%
  full_join(., n_peaks, by = c('ASV', 'Farm')) %>% # join number of peaks data
  mutate(# quality of assigments
    Quality = as.numeric(as.character(cut(
      diff,
      breaks = c(0, 100, 200, 1000), # difference where ASV reads peak in relation to TFS 
      labels = c(1, 2, 3)
    ))),
    Quality = if_else(n_peaks > 2, 4, Quality)) # double peaks)


# # Quality filteriing ----
top_qual_ASV <- 
  peak_dat %>%
  select(ASV, Farm, Quality) %>%
  spread(Farm, Quality) %>%
  drop_na() %>% 
  filter(Ext<3 & Int <3) %>% 
  select(ASV) %>% 
  as.list()

## join peak_dat with splines and ASV reads abundance data --------------
long_top_ASV <- 
    full_join(bact_top_long, peak_dat, by = c('ASV', 'Farm')) %>% 
  dplyr::filter(ASV %in% top_qual_ASV$ASV)
  
# EG assigment ----------
EG_groups <-
  peak_dat %>%
  dplyr::filter(ASV %in% top_qual_ASV$ASV) %>%
  group_by(ASV) %>%
  summarise(mean_peak = mean(Peak)) %>%
  mutate(EG = cut(
    mean_peak,
    breaks = c(0, quantile(data_all$TFS, probs = c(1 / 5, 2 / 5, 3 / 5, 4 / 5)), max(data_all$TFS)),
    # breaks = c(0, 140.8, 140*2, 140.8*3, 140.8*4, 140.8 * 5),
    labels = c("I", "II", "III", "IV", "V")
  )) %>%
  dplyr::select(ASV, EG) %>%
  write_csv('outputs/EG_groups_MED.csv')

# print the table with EG groups

kable(EG_groups,
      caption = "EcoGroup assigment for the 2,000 most abundant bacterial ASVs based on quantile spline regressions against ES")  

# plot spline examples of EG I - V ----------
long_top_ASV_EG <-
  full_join(long_top_ASV, EG_groups, by = "ASV") %>%
  filter(n_peaks == 1 & Quality == 1) %>%
  group_by(EG, ASV) %>%
  summarise(diff = first(diff)) %>% top_n(n = 2)

eg_example_plot <- 
  long_top_ASV %>% 
  filter(ASV %in% long_top_ASV_EG$ASV) %>% 
  ggplot(aes(x = TFS, y = abund, color = Farm)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'Total free sulfides', y = 'Number of reads') +
  theme_javier() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  geom_vline(aes(xintercept = Peak, color = Farm)) +
  geom_vline(aes(xintercept = Peak, color = Farm)) +
  facet_wrap(~ ASV, scales = 'free') 

print(eg_example_plot)
ggsave(eg_example_plot, 
       filename = 'figures/eg_example_plot.tiff',
       width = 8,
       height = 6,
       dpi = 300,
       device = 'tiff',
       compression = 'lzw')



# proportion of AVS by EG groups-----
EG_groups %>% 
  group_by(EG) %>% 
  summarise(n = n()/nrow(.)*100)


# Summary plots of eg assigments ---------------
pp1 <- 
  ggplot(EG_groups, aes(fct_rev(EG))) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), color = 1) + 
  scale_y_continuous(labels=scales::percent) +
  labs(y  = "Relative frequencies", x = "Eco_Group") +
  coord_flip()

pp2 <- 
  ggplot(peak_dat, aes(fct_rev(factor(Quality)))) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), color = 1) + 
  scale_y_continuous(labels=scales::percent) +
  labs(y  = "Relative frequencies", x = "Spline Quality Score") +
  coord_flip()

ggarrange(pp1,pp2)

ggsave(ggarrange(pp1,pp2), 
       filename = 'figures/EG_assigment_summary.tiff',
       device = 'tiff',
       compression = 'lzw',
       width = 8,
       height = 4,
       dpi = 300)


# Calculate the bacterial AMBI------------------

# read AVS read abundance and eco-groups data-----------------------------------
data <-  
  data_all %>%
  select(-c(2:12))

# bact_data_pa <- 
#   data %>% 
#   mutate_if(is.numeric,~if_else(.>0,1,0))


# EG proportion calculations --------------
eg_prop <- 
data %>% 
  gather(ASV, abund, -SampleID) %>%
  inner_join(EG_groups, by = 'ASV') %>%
  group_by(SampleID) %>%
  mutate(N = sum(abund)) %>% 
  group_by(SampleID, EG) %>%
  mutate(n = sum(abund) / N) %>%
  summarise(n = first(n) * 100) %>%
  spread(EG, n)


# Combine all indices with metadata----------------
env <- 
  data_all %>% 
  select(SampleID:TOC) %>% 
  write_csv('data/env_bact.csv')

EG_prop_all <- left_join(eg_prop, env, by = 'SampleID')

# Calculate optimal weights for AMBI using linear regression on the train dataset-------
AMBI_weights <- lm(TFS ~  0 + I + II + III + IV + V, data = EG_prop_all)

indices_all <-
  EG_prop_all %>%
  mutate(
    bMBI = (
      AMBI_weights$coefficients[1] * I +
        AMBI_weights$coefficients[2] * II + 
        AMBI_weights$coefficients[3] * III + 
        AMBI_weights$coefficients[4] * IV + 
        AMBI_weights$coefficients[5] * V
    )/100
  ) %>% 
  ungroup() %>% 
  mutate(Distance = fct_relevel(Distance, "0m", "50m", "100m", "250m", "500m", "1000m","2000m", "Control" ))


indices_all %>%
  select(SampleID, bMBI) %>%
  rename(AMBI_bact = bMBI) %>%
  write_csv('outputs/bMBI_data.csv')

# Plot euk-AMBI vs macrofaunal AMBI-----------
plot4 <-
  ggplot(indices_all, aes(TFS, bMBI)) +
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

print(plot4)

ggsave(plot4, 
       filename = 'figures/bMBI_vs_TFS.tiff',
       device = 'tiff',
       compression = 'lzw',
       width = 6,
       height = 4,
       dpi = 300)