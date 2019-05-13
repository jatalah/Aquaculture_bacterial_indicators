# load libraries------------
library(tidyverse)
library(splines)
library(quantreg)
library(knitr)
source('theme_javier.R')
theme_set(theme_javier())

#read salmon farm foram DNA only data --------------
dna <- 
  read_csv('BacteriaMed_all.csv',col_types = cols())

# select 150 most abundant forams OTUs------------
abund_OTU <- 
  dna %>% 
  gather(OTU, abund, OTU1:OTU12955) %>% 
  group_by(OTU) %>% 
  summarise_at(vars(abund), funs(sum), na.rm =T) %>%
  top_n(200) %>% 
  dplyr::select(OTU) %>% 
  as.list()

dna_top <- 
  dna %>% 
  select(Station, Distance, Installation, TFS, abund_OTU$OTU)

# convert into long format---------------
dna_top_long <- 
  dna_top %>% 
  gather(OTU,abund, abund_OTU$OTU)

# create plot with quantile regression and splines with 3 df--------------
p2 <- 
  ggplot(dna_top_long, aes(x = Distance, y = abund, color = Installation)) +
  geom_point(size = 3, alpha = .3) +
  labs(x = 'Distancia (m)', y = 'Number of reads') +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~OTU, scales = 'free')


# get splines data from plot------------------------
gb <- ggplot_build(p2)

# obtaine ES values whre abundance peak by year-------
gb1 <- 
  gb$data[[2]] %>% 
  data.frame() %>% 
  group_by(group,PANEL) %>% 
  nest() %>% 
  mutate(peak = map(data, ~.$x[which.max(.$y)]))

peak_dat <- 
  gb1$peak %>% 
  flatten_dbl() %>% 
  data.frame() %>% 
  rename(Peak = '.') %>% 
  bind_cols(dna_top_long %>% distinct(OTU,Installation)) %>% 
  spread(Installation,Peak) %>% 
  mutate(diff = abs(Farm2 - Farm1)) %>% 
  filter(diff<100) %>% 
  gather(Installation, Peak, Farm2:Farm1) %>% 
  select(-diff)


long_top_OTU <- 
  inner_join(dna_top_long, peak_dat, by = c('OTU', 'Installation'))
  
## plot selected OTUs with peak vertical lines--------------
plot_distance <-
  ggplot(long_top_OTU, aes(x = Distance, y = abund, color = Installation)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'Distance to the farms (m)', y = 'Number of reads') +
  theme_javier() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ OTU, scales = 'free') +
  geom_vline(aes(xintercept = Peak, color = Installation))

ggsave(plot_distance, 
       filename = 'Distance_MED.tiff',
       width = 15,
       height = 10,
       device = 'tiff',
       compression = 'lzw')

## Sulphide gradient ---------------
# create plot with quantile regression and splines with 3 df--------------
p3 <- 
  ggplot(dna_top_long, aes(x = TFS, y = abund, color = Installation)) +
  geom_point(size = 3, alpha = .3) +
  labs(x = 'Distancia (m)', y = 'Number of reads') +
  theme_bw() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~OTU, scales = 'free')


# get splines data from plot------------------------
gb_tfs <- ggplot_build(p3)

# obtaine ES values whre abundance peak by year-------
gb_tfs1 <- 
  gb_tfs$data[[2]] %>% 
  data.frame() %>% 
  group_by(group,PANEL) %>% 
  nest() %>% 
  mutate(peak = map(data, ~.$x[which.max(.$y)]))

peak_dat_tfs <- 
  gb_tfs1$peak %>% 
  flatten_dbl() %>% 
  data.frame() %>% 
  rename(Peak = '.') %>% 
  bind_cols(dna_top_long %>% distinct(OTU,Installation)) %>% 
  spread(Installation,Peak) %>% 
  mutate(diff = abs(Farm2 - Farm1)) %>% 
  filter(diff<100) %>% 
  gather(Installation, Peak, Farm2:Farm1) %>% 
  select(-diff)


long_top_OTU_tfs <- 
  inner_join(dna_top_long, peak_dat_tfs, by = c('OTU', 'Installation'))

## plot selected OTUs with peak vertical lines--------------
plot_tfs <- 
  ggplot(long_top_OTU_tfs, aes(x = TFS, y = abund, color = Installation)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'Total free sulfides', y = 'Number of reads') +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap( ~ OTU, scales = 'free') +
  geom_vline(aes(xintercept = Peak, color = Installation))

ggsave(plot_tfs, 
       filename = 'TFS_MED.tiff',
       width = 15,
       height = 10,
       device = 'tiff',
       compression = 'lzw')

print(plot_tfs)

# EG assigment ----------
EG_groups_MED <-
  peak_dat_tfs %>%
  mutate(EG = cut(
    Peak,
    breaks= c(0,quantile(dna$TFS,probs = c( 1/4,2/4,3/4)),max(dna$TFS)),
    labels = c("I", "II", "III","IV")
  )) %>% 
  select(OTU,EG) %>% 
  write_csv(.,'EG_groups_MED.csv')

# print the table with EG groups-------------

kable(breaks)

kable(EG_groups_MED, 
      caption = "EcoGroup assigment for the 100 most abundant foraminifera OTUs based on quantile spline regressions against ES")  
