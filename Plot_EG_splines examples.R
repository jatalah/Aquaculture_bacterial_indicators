library(tidyverse)
library(quantreg)

long_top_ASV <- read_csv('data/long_top_ASV.csv')

# ASV selected by Eva for each EG I - V in order----
spline_examples <- c("ASV_1115","ASV_475", "ASV_321", "ASV_348", "ASV_169")

eg_example_plot <-
  long_top_ASV %>%
  dplyr::filter(ASV %in% spline_examples) %>%
  left_join(EG_groups, by = "ASV") %>%
  mutate(Farm = fct_recode(Farm, `Farm 1` = "Ext", `Farm 2` = "Int")) %>% 
  ggplot(aes(x = TFS, y = abund, color = Farm)) +
  geom_point(size = 2, alpha = .3) +
  labs(x = 'Acide Volatile Sulphides'~"("*mu*M*")", y = 'Number of reads', parse = T) +
  theme_javier() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  geom_vline(aes(xintercept = Peak, color = Farm)) +
  geom_vline(aes(xintercept = Peak, color = Farm)) +
  facet_wrap(~ ASV, scales = 'free', ncol = 5) +
  scale_color_manual(values = c("#A95AA1", "#F5793A"))

print(eg_example_plot)

# save plot----------------
ggsave(
  eg_example_plot,
  filename = 'figures/eg_example_plot.tiff',
  width = 10.5,
  height = 2.45,
  dpi = 300,
  device = 'tiff',
  compression = 'lzw'
)

ggsave(
  eg_example_plot,
  filename = 'figures/eg_example_plot.svg',
  width = 10.5,
  height = 2.45,
  dpi = 300
)


# plot spline examples of EG I - V ----------
# long_top_ASV_EG <-
#   full_join(long_top_ASV, EG_groups, by = "ASV") %>%
#   dplyr::filter(n_peaks == 1 & Quality > 1) %>%
#   group_by(EG, ASV) %>%
#   summarise(diff = first(diff)) %>%
#   top_n(n = 5)
