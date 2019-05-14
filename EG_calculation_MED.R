# load libraries------------
library(tidyverse)
library(splines)
library(quantreg)
library(knitr)
library(pracma)# for finding peaks
library(ggpubr)
library(vegan)
library(ggpmisc)

source('theme_javier.R')
theme_set(theme_javier())

#read data --------------
dna <- read_csv('data/Bacteria16SMed_Table.csv',col_types = cols())

# select 500 most abundant ASV------------
abund_ASV <-  
  dna %>% 
  gather(ASV, abund, ASV_001:ncol(.)) %>% 
  group_by(ASV) %>% 
  summarise_at(vars(abund), funs(sum), na.rm =T) %>%
  top_n(2000) %>% 
  dplyr::select(ASV) %>% 
  as.list()

dna_top <- 
  dna %>% 
  select(Station, Distance, Installation, TFS, abund_ASV$ASV)

# convert into long format
dna_top_long <- 
  dna_top %>% 
  gather(ASV,abund, abund_ASV$ASV)

# create plot with quantile regression and splines with 3 df--------------
p2 <- 
  ggplot(dna_top_long, aes(x = TFS, y = abund, color = Installation)) +
  geom_point(size = 3, alpha = .3) +
  labs(x = 'Total free sulphides', y = 'Number of reads') +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ASV, scales = 'free')


# get splines data from plot------------------------
gb <- ggplot_build(p2)

# obtaine ES values where abundance peak by installation-------
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
  distinct(Installation,ASV) %>% 
  arrange(ASV)


# get the number of peaks per AVS and year ------
n_peaks <- 
  gb1 %>% 
  select(n_peaks) %>% 
  unnest() %>% 
  bind_cols(plot_id)


# merge difference peaks by Installation  and the number of peaks-----
peak_dat <-
  gb1 %>%
  select(peak) %>%
  unnest() %>%
  bind_cols(plot_id) %>%
  spread(Installation, peak) %>%
  mutate(diff = abs(Ext - Int)) %>%
  gather(Installation, Peak, Ext:Int) %>%
  full_join(., n_peaks, by = c('ASV', 'Installation')) %>% # join number of peaks data
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
  select(ASV, Installation, Quality) %>%
  spread(Installation, Quality) %>%
  drop_na() %>% 
  filter(Ext<3 & Int <3) %>% 
  select(ASV) %>% 
  as.list()

## join peak_dat with splines and ASV reads abundance data --------------
long_top_ASV <- 
  full_join(dna_top_long, peak_dat, by = c('ASV', 'Installation')) %>% 
  dplyr::filter(ASV %in% top_qual_ASV$ASV)
  
## plot selected ASVs with peak vertical lines--------------
plot_TFS <-
  ggplot(long_top_ASV, aes(x = TFS, y = abund, color = Installation)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'Total free sulfides', y = 'Number of reads') +
  theme_javier() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  geom_vline(aes(xintercept = Peak, color = Installation)) +
  facet_wrap(~ ASV, scales = 'free') 

ggsave(plot_TFS, 
       filename = 'figures/TFS_splines.tiff',
       width = 25,
       height = 20,
       dpi = 90,
       device = 'tiff',
       compression = 'lzw')

# EG assigment ----------
EG_groups <- 
  peak_dat %>%
  dplyr::filter(ASV %in% top_qual_ASV$ASV) %>%  
  group_by(ASV) %>% 
  summarise(mean_peak = mean(Peak)) %>% 
  mutate(EG = cut(
    mean_peak,
    breaks= c(0, quantile(dna$TFS,probs = c( 1/5,2/5,3/5, 4/5)),max(dna$TFS)),
    labels = c("I", "II", "III","IV", "V")
  )) %>% 
  dplyr::select(ASV,EG) %>% 
  write_csv('outputs/EG_groups_MED.csv')

# print the table with EG groups-------------

kable(EG_groups, 
      caption = "EcoGroup assigment for the 100 most abundant foraminifera ASVs based on quantile spline regressions against ES")  


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


ggsave(ggarrange(pp1,pp2), 
       filename = 'figures/EG_assigment_summary.tiff',
       device = 'tiff',
       compression = 'lzw',
       width = 8,
       height = 4,
       dpi = 300)


# Calculate the bacterial AMBI------------------

# read AVS read abundance and eco-groups data-----------------------------------
env <-  
  dna %>%
  select(Station:TOC)

data <- 
  dna %>%
  select(Station, ASV_001:ncol(.))

# EG proportion calculations --------------
EG_indices <-
  data %>%
  gather(ASV, abund, -Station) %>%
  left_join(EG_groups, by = 'ASV') %>%
  group_by(Station) %>%
  mutate(
    N = sum(abund),
    H = diversity(abund),
    S = specnumber(abund),
    J = H / log(S),
    UA = sum(is.na(EG) & abund > 0) / S * 100
  ) %>%
  drop_na(EG) %>%
  mutate(N1 = sum(abund),
         S1 = specnumber(abund)) %>%
  group_by(Station, EG) %>%
  mutate(prop_n = sum(abund) / N1,
         prop_s = specnumber(abund) / S1) %>%
  summarise(
    n = mean(prop_n, na.rm = T) * 100,
    s = mean(prop_s, na.rm = T) * 100,
    H = mean(H),
    S = mean(S),
    J = mean(J),
    N = mean(N),
    N1 = mean(N1),
    S1 = mean(S1),
    UA = mean(UA)
  ) %>%
  gather(key, value, n, s) %>%
  unite(EG_comb, EG, key) %>%
  spread(EG_comb, value) 


# Combine all indices with metadata----------------
EG_prop_all <- left_join(EG_indices, env, by = 'Station')

# Calculate optimal weights for AMBI using linear regression on the train dataset-------
AMBI_weights <- lm(TFS ~  0 + I_n + II_n + III_n + IV_n + V_n, data = EG_prop_all)

indices_all <-
  EG_prop_all %>%
  mutate(
    AMBI1 = (
      1.5 * II_n + 
        3 * III_n + 
        4.5 * IV_n + 
        6 * V_n
    ),
    AMBI2 = (
      AMBI_weights$coefficients[1] * I_n +
        AMBI_weights$coefficients[2] * II_n + 
        AMBI_weights$coefficients[3] * III_n + 
        AMBI_weights$coefficients[4] * IV_n + 
        AMBI_weights$coefficients[5] * V_n
    )/100
  ) %>% 
  ungroup() %>% 
  mutate(Distance = fct_relevel(Distance, "0m", "50m", "100m", "250m", "500m", "1000m","2000m", "Control" ))

# Plot euk-AMBI vs macrofaunal AMBI-----------
plot4 <- 
  ggplot(indices_all, aes(TFS, AMBI2)) +
  geom_point(size = 4, aes(color  = Distance, shape = Installation)) +
  geom_smooth(method = "lm") +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = ..rr.label..),
    parse = TRUE,
    size = 4
  ) +
  labs(y = 'Bacterial-AMBI', x =  'Total free sulphide (TFS)')

ggsave(plot4, 
       filename = 'figures/bct_AMBI_vs_TFS.tiff',
       device = 'tiff',
       compression = 'lzw',
       width = 6,
       height = 4,
       dpi = 300)
