# load libraries------------
library(tidyverse)
library(splines)
library(quantreg)
library(knitr)
library(pracma)# for finding peaks
library(ggpubr)
library(vegan)
library(phyloseq)
library(broom)
library(ggfortify)

source('theme_javier.R')
theme_set(theme_javier())

bact_top_long <- read_csv('data/bact_top_long.csv')

# create plot with quantile regression and splines with 3 df--------------
p2 <-
  ggplot(bact_top_long, aes(x = TFS, y = abund, color = Farm)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ ASV, scales = 'free')


# get splines data from plot------------------------
gb <- ggplot_build(p2)

# obtaine ES values where abundance peak by Farm-------
gb1 <-
  gb$data[[2]] %>%
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0,
        minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )


plot_id <-
  gb$plot[[1]] %>%
  distinct(Farm, ASV) %>%
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
  ungroup() %>%
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
      breaks = c(0, 25, 50, 100),
      # difference where ASV reads peak in relation to TFS
      labels = c(1, 2, 3)
    ))),
    Quality = if_else(n_peaks > 2, 4, Quality)) # double peaks)

## Quality filteriing ----
top_qual_ASV <-
  peak_dat %>%
  select(ASV, Farm, Quality) %>%
  spread(Farm, Quality) %>%
  drop_na() %>%
  filter(Ext < 3 & Int < 3) %>%
  select(ASV) %>%
  as.list()

## join peak_dat with splines and ASV reads abundance data --------------
long_top_ASV <-
  full_join(bact_top_long, peak_dat, by = c('ASV', 'Farm')) %>%
  dplyr::filter(ASV %in% top_qual_ASV$ASV) %>% 
  write_csv('data/long_top_ASV.csv')

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


# print the table with EG groups-------------
kable(EG_groups,
      caption = "EcoGroup assigment for the 2,000 most abundant bacterial ASVs based on quantile spline regressions against ES")

# proportion of AVS by EG groups-----
EG_groups %>%
  group_by(EG) %>%
  summarise(n = n() / nrow(.) * 100)


# Summary plots of eg assigments ---------------
pp1 <-
  ggplot(EG_groups, aes(fct_rev(EG))) +
  geom_bar(aes(y = (..count..) / sum(..count..)), color = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(y  = "Relative frequencies", x = "Eco_Group") +
  coord_flip()

pp2 <-
  ggplot(peak_dat, aes(fct_rev(factor(Quality)))) +
  geom_bar(aes(y = (..count..) / sum(..count..)), color = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(y  = "Relative frequencies", x = "Spline Quality Score") +
  coord_flip()

ggarrange(pp1, pp2)

ggsave(
  ggarrange(pp1, pp2),
  filename = 'figures/EG_assigment_summary.tiff',
  device = 'tiff',
  compression = 'lzw',
  width = 8,
  height = 4,
  dpi = 300
)