

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
path <- '../data/tif/index'
dirs <- dir_ls(path, type = 'directory')


# Baseline ----------------------------------------------------------------

fles.bsln <- grep('baseline', dirs, value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  map(., dir_ls) %>% 
  unlist() %>% 
  as.character()
fles.bsln

