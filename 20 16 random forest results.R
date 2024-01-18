

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue, gtools, stringr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

dirs <- dir_ls('../data/tif/models_dataverse')
dirs <- as.character(dirs)

## Colombia ------------------------------
fles.col <- grep('col', dirs, value = T) %>% dir_ls(., regexp = '.tif$') %>% as.character()
rfrs.col <- grep('RF', fles.col, value = T) %>% map(., rast)
rfrs.col[[3]] <- terra::crop(rfrs.col[[3]], ext(rfrs.col[[1]])) %>% terra::mask(., rfrs.col[[1]])
rfrs.col <- reduce(rfrs.col, c)

impr.col <- grep('IG', fles.col, value = T) %>% map(., rast)
impr.col <- reduce(impr.col, c)

### Labels Colombia 
lbls.col <- tibble(
  value = 1:7, 
  class = c('Sin idoneidad', 'Sin idoneidad', 'Caliente - Húmedo', 'Templado - Seco', 'Muy caliente - Húmedo', 'Limitaciones', 'Aptidud incierta')
)

## Ecuador ------------------------------
fles.ecu <- grep('ecu', dirs, value = T) %>% dir_ls(., regexp = '.tif$') %>% as.character()
rfrs.ecu <- grep('RF', fles.ecu, value = T) %>% rast()

### Labels Ecuador 
lbls.ecu <- tibble(
  value = 1:8,
  class = c('Sin idoneidad', 'Sin idoneidad', 'Templado - Seco', 'Cálido - Húmedo', 'Frío - Seco', 'Cálido - Seco', 'Limitaciones', 'Aptitud incierta')
)

## Perú ---------------------------------
fles.per <- grep('per', dirs, value = T) %>% dir_ls(., regexp = '.tif$') %>% as.character()
rfrs.per <- grep('RF', fles.ecu, value = T) %>% rast()

### Labels Perú
lbls.per <- tibble(
  value = 1:9, 
  class = c('Sin idoneidad', 'Sin idoneidad', 'Templado - Seco', 'Cálido - Húmedo', 'Muy cálido - Seco', 'Frío - Seco', 'Cálido - Seco', 'Limitaciones', 'Aptitud incierta')
)

# To write the labels -----------------------------------------------------
write.csv(lbls.col, '../data/tif/models_dataverse/labels_aez-col.csv', row.names = FALSE)
write.csv(lbls.ecu, '../data/tif/models_dataverse/labels_aez-ecu.csv', row.names = FALSE)
write.csv(lbls.per, '../data/tif/models_dataverse/labels_aez-per.csv', row.names = FALSE)

# Labels impact gradient 
lbls.imp <- tibble(
  value = c(0, 1, 3, 4, 5),
  impact = c('Sin idoneidad', 'Adaptaicón incremental', 'Transformación', 'Oportunidades', 'Resiliencia sistémica')
)

write.csv(lbls.imp, '../data/tif/models_dataverse/labels_imp.csv', row.names = FALSE)


