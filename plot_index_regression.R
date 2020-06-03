library(tidyverse)
library(ggpmisc)
indices_all <- read_csv('data/indices_all.csv')
source('theme_javier.R')
theme_set(theme_javier())

# Plot euk-AMBI vs macrofaunal AMBI-----------
lm_index <- lm(bMBI_keely ~ TFS, indices_all)
tidy(summary(lm_index))
tidy(summary(lm_index))
autoplot(lm_index) + theme_bw()

plot4 <-
  ggplot(indices_all,
         aes(TFS, bMBI_keely)) +
  geom_point(size = 4, aes(color  = Distance, shape = Farm),alpha = .8) +
  stat_smooth(method = "lm",
              formula = y ~ sqrt(x),
              alpha = .2) +
  stat_poly_eq(
    formula = y ~ sqrt(x),
    aes(label =  stat(rr.label)),
    parse = TRUE,
    size = 4
  ) +
  labs(y = 'Bacterial Biotic Index',
       x = 'Acide Volatile Sulphides' ~ "(" * mu * M * ")",
       parse = T) +
  scale_color_manual(
    values = c(
      "lightcoral", # 0 m 
      "firebrick2", # 50 
      "goldenrod", # 100 m
      "gold4", # 250
      "saddlebrown",# 500 m
      "lightskyblue1", # 1000 m 
      "royalblue1", # 2000 m
      "navy" # Control
    )
  )

print(plot4)

ggsave(
  plot4,
  filename = 'figures/bMBI_vs_TFS.tiff',
  device = 'tiff',
  compression = 'lzw',
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  plot4,
  filename = 'figures/bMBI_vs_TFS.svg',
  device = 'svg',
  width = 6,
  height = 4,
  dpi = 300
)
