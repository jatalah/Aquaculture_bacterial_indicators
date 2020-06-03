theme_javier <- function(base_size = 12, base_family = "")
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    )
}

