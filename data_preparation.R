# load libraries------------
library(tidyverse)

# EXT is Farm 1; INT is Farm 2.
# El Código de colores para FARM 1 es #A95AA1 y para Farm 2  #F5793A

#read data --------------
data_all <- read_csv('data/bacteria.csv', col_types = cols())

# data transformatiuons if required ----
data_pa <- data_all %>% mutate_if(is.numeric,  ~ if_else(. > 0, 1, 0))
data_sqrt <- data_all %>% mutate_if(is.numeric,  ~ sqrt(.))

env <-
  data_all %>%
  select(SampleID:TOC) %>%
  write_csv('data/env_bact.csv')


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
  left_join(env, ., by = 'SampleID')

# select 2000 most abundant ASV------------
abund_ASV <-
  bact_rare_df %>%
  gather(ASV, abund, ASV_001:ncol(.)) %>%
  group_by(Farm, ASV) %>%
  summarise_at(vars(abund), list(sum), na.rm = T) %>%
  top_n(n = 1000, wt = abund) %>%
  ungroup() %>%
  # complete(ASV, nesting(Farm)) %>% # to remove ASV not present in both farms
  # group_by(ASV) %>%
  # filter(!is.na(sum(abund))) %>%
  # ungroup() %>%
  select(ASV) %>%
  as.list()

bact_top <-
  bact_rare_df %>%
  select(Station, Distance, Farm, TFS, abund_ASV$ASV)

# convert into long format
bact_top_long <-
  bact_top %>%
  gather(ASV, abund, abund_ASV$ASV) %>% 
  write_csv('data/bact_top_long.csv')
