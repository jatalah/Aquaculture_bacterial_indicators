library(tidyverse)
library(maps)
library(mapdata)
library(ggrepel)
library(sf)


coords <- read_csv('data/Coordinates_Spain.csv')

ggplot(coords,aes(x, y, color = Location, label = Station)) +
  geom_point() +
  theme_bw() +
  geom_label_repel() +
  labs(x = 'Lon' , y = 'Lat')


# calculate distance between stations----------
# Instalation int---
coord_int <-
  coords %>%
  filter(Location=='Int') %>% 
  st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) 

centroid_int<- 
  coord_int %>%
  rownames_to_column() %>% 
  filter(rowname==12) %>% 
  st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) 


dist_int <- 
  coord_int %>% 
  mutate(distance = st_distance(coord_int, centroid_int)) %>% 
  as.data.frame()


# Instalation ext----
coord_ext <-
  coords %>%
  filter(Location=='Ext') %>% 
  st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) 

centroid_ext<- 
  coord_ext %>%
  rownames_to_column() %>% 
  filter(rowname==12) %>% 
  st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) 

dist_ext <- 
  coord_ext %>% 
  mutate(distance = st_distance(coord_ext, centroid_ext)) %>% 
  as.data.frame() 

x <- bind_rows(dist_int,dist_ext)
