# load libraries------------
library(tidyverse)
library(splines)
library(ggpubr)
library(vegan)
library(ggord)
library(ggpmisc)
library(readxl)
library(knitr)
library(emmeans)
library(pairwiseAdonis)
library(benthos)

source('theme_javier.R')
theme_set(theme_javier())

# read data -----
poli <- 
  read_excel('data/poliquetos.xlsx', 2) %>% 
  mutate(Distance = fct_relevel(Distance, "0m" , "50m", "250m", "2000m"))

factors <- 
  poli %>% 
select(SampleID:Rep)

env <- 
  read_excel('data/poliquetos.xlsx', 1) %>% 
  select(SampleID, TN, TC, MO, TFS)


poli_env <- left_join(poli, env, by = 'SampleID')


poli_long <- 
  poli %>% 
  gather(taxa, abu, Capitellidae:Nephtyidae) 

top_taxa <- 
  poli_long %>% 
  group_by(taxa) %>% 
  summarise_at(vars(abu), funs(sum), na.rm =T) %>%
  top_n(12) %>% 
  dplyr::select(taxa) %>% 
  as.list()

poli_long %>%
  filter(taxa %in% top_taxa$taxa) %>% 
  ggplot() +
  geom_boxplot(aes(Distance, abu, color = Farm)) +
  facet_wrap( ~ taxa, scales = 'free', nrow = 3)

# PERMANOVA---------
permanova <-
  adonis(
    sqrt(sqrt(poli[,-c(1:5)])) ~ Distance * Farm,
    method = 'bray',
    permutations = 9999,
    data = poli
  )

permanova


x <- 
  adonis(
    sqrt(sqrt(poli_env[,c(6:26)]))  ~ TC + MO + TFS + TN,
    method = 'bray',
    permutations = 999,
data = poli_env,
by = 'margin'
  )

adonis2(
  sqrt(sqrt(poli_env[,c(6:26)]))  ~ TFS + TC + MO +  TN,
  method = 'bray',
  # permutations = 999,
  data = poli_env,
  by = "margin"
)

dist_lm <- 
  adonis2(
  sqrt(sqrt(poli_env[,c(6:26)]))  ~ TFS + TC + MO +  TN,
  method = 'bray',
  permutations = 999,
  data = poli_env,
  by = "margin"
)

AICc.PERMANOVA(adonis.model = x)


AICc.PERMANOVA(x)

# 08 PERMANOVA pairwise comparisons-------------
pairwise.adonis(
  sqrt(sqrt(poli[,-c(1:5)])),
  poli$Distance)

# 04 MDS ordination -------------
mds <- metaMDS( sqrt(sqrt(poli[, -c(1:5)])) , distance = 'bray')

## MDS biplot 
ggord(
  mds,
  grp_in = poli$Distance,
  poly = F,
  alpha = 1,
  ellipse = T,
  arrow = 0,
  repel = T,
  text = .01,
  vec_ext = .7
) +
  theme_javier()

# SIMPER ---------------
simp_distance <- simper(sqrt(sqrt(poli[, -c(1:5)])), poli$Distance, permutations = 999)
simp_farm <- simper(sqrt(sqrt(poli[, -c(1:5)])), poli$Farm, permutations = 999)


simper_dist_table <- 
  bind_rows("0m vs. 50m" = simp_distance$`0m_50m` %>% data.frame(),
            "0m vs 250m" = simp_distance$`0m_250m` %>% data.frame(),
            "0m vs 2000m" = simp_distance$`0m_2000m` %>% data.frame(),
          "50m vs 250m" = simp_distance$`50m_250m` %>% data.frame(),
          "50m vs. 2000m" = simp_distance$`50m_2000m` %>% data.frame(),
          "250m vs 2000m" = simp_distance$`250m_2000m` %>% data.frame(),
            .id = "Groups") %>% 
  dplyr::select(-overall, -ord) %>% 
  dplyr::filter(cusum<.7) %>% 
  rename(Family = "species")
  
kable(simper_dist_table, digits = 2)

# PERMDIST----
dis <- vegdist(sqrt(sqrt(poli[, -c(1:5)])))
beta_disp <-
  betadisper(
    dis, poli$Distance
  )
beta_disp
anova(beta_disp)
permutest(beta_disp, pairwise = TRUE, permutations = 999)
plot(beta_disp)
boxplot(beta_disp)

# Univariate diversity index -------------------
indices <- 
  poli %>% 
  transmute(N = rowSums(poli[, -c(1:5)]),
            H = diversity(poli[, -c(1:5)]),
            S = specnumber(poli[, -c(1:5)]),
            J = H/log(S)) %>% 
  bind_cols(poli[, c(1:5)]) %>% 
  write_csv('outputs/univariate_poli_indices.csv')



boxplots <- 
  indices %>% 
  gather(index, value, c("N", "S", "H", "J")) %>% 
  mutate(index = fct_relevel(index, c("N", "S", "J"))) %>% 
  ggplot(., aes(x = Distance, y = value, fill = Farm)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~index,scales = 'free') +
  theme_bw()

boxplots


indices_dat <-
  indices %>%
  mutate(N = log(N +1)) %>% 
  gather(index, value, c("N", "S", "H", "J")) %>%
  group_by(index) %>%
  nest() %>%
  mutate(
    anova = map(.x = data, ~ aov(value ~ Distance * Farm, data = .x)),
    anova_tab = map(anova, ~car::Anova(.)),
    anova_table = map(anova_tab, broom::tidy),
    residuals = map(anova, broom::augment))

anova_table <- 
  indices_dat %>% 
  select(index,anova_table) %>% 
  unnest()

kable(anova_table,caption = 'Univaraite ANOVAs')

## model validation----
res <- 
  indices_dat %>% 
  select(index, residuals) %>% 
  unnest(residuals, .drop = TRUE)

# Plot of fitted vs. residual values by index to check the assumptions heteroscedasticity or any other pattern,
#  e.g. non-linearity or missing covariates
ggplot(res) +
  geom_point(aes(x = .fitted, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2) 

# boxplot by type
ggplot(res) +
  geom_boxplot(aes(x = Distance, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# qqplot of the normalised residuals to check the assumption of normality------------
ggplot(res) +
  stat_qq(aes(sample = .std.resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_abline(
    intercept =  0,
    slope = 1,
    lty = 2,
    col = 2
  )

# calculate AMBI------------
poli_long %>% 
  mutate(COMPLIANT = is_accepted(taxon = taxa)) %>% 
  group_by(taxa) %>% 
  distinct(COMPLIANT) %>% 
  filter(COMPLIANT==FALSE) # two non-compliant taxa

AMBI_poli <- 
  poli_long %>%
  group_by(SampleID) %>% 
  nest() %>% 
  mutate(AMBI = map(data, ~ ambi(taxon = .x$taxa, count = .x$abu))) %>% 
  select(SampleID,AMBI) %>% 
  unnest() %>% 
  right_join(indices, by = "SampleID")

ggplot(AMBI_poli, aes(x = Distance, y = AMBI, fill = Farm)) +
  geom_boxplot(alpha = .8)

bact_ambi <- read_csv('data/bact_AMBI.csv')
bact_poli_ambi <- left_join(AMBI_poli, bact_ambi, by = "SampleID")


bact_poli_ambi %>% 
  filter(SampleID !="Int_2000m_1_c") %>% 
  ggplot(aes(x = AMBI, y = AMBI_bact)) +
  geom_point(alpha = .8, size = 5, aes(color = Farm)) +
  geom_smooth(method = 'lm') +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = ..rr.label..),
    parse = TRUE,
    size = 4
  )

bact_poli_ambi %>% 
  filter(SampleID !="Int_2000m_1_c") %>% # seems to be an outlier with 82 Capetilidae at 200 m

ggplot(aes(x = AMBI, y = AMBI_bact)) +
  geom_point(alpha = .8, size = 5, aes(color = Farm)) +
  geom_smooth(method = 'lm') +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = ..rr.label..),
    parse = TRUE,
    size = 4
  )


ggplot(bact_poli_ambi, aes(x = AMBI, y = AMBI_bact, label = SampleID)) +
  geom_text() 

AMBI_env <- left_join(AMBI_poli, env, by = "SampleID")

# AMBI enviroment relationships ------------
AMBI_env %>% 
  gather(key, value, TN:TFS) %>% 
ggplot(aes(x = value, y = AMBI)) +
  geom_point(alpha = .8, aes(color = Farm), size = 5) +
  geom_smooth() +
  facet_wrap(~key, scales = 'free')