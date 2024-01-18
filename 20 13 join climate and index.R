

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Raster data
clma <- terra::rast('../data/tif/results/climate_stack_col-ecu-per.tif')
indx <- terra::rast('../data/tif/results/indx_bsl-ftr_chirts-ssp370.tif')

# Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- wrld[wrld$sov_a3 %in% c('COL', 'ECU', 'PER'),]
zone <- vect(zone)

# Join between only one stack ---------------------------------------------

# Create the masks
mask.clma <- (clma[[1]] * 0) + 1
mask.clma <- setNames(object = mask.clma, nm = 'mask.clma')
mask.indx <- (indx[[1]] * 0) + 1
mask.indx <- setNames(object = mask.indx, nm = 'mask.indx')

# Crop 
indx <- terra::crop(indx, ext(mask.clma))


# Raster to table ---------------------------------------------------------

clma.tble <- as_tibble(terra::as.data.frame(clma, xy = T))
nrow(clma.tble)

# Index to table
indx.tble <- as_tibble(terra::as.data.frame(indx, xy = T))
indx.tble <- drop_na(indx.tble)
indx <- terra::rast(indx.tble, type = 'xyz')
nrow(indx.tble)

# Create the mask for the index
mask.indx <- indx[[1]] * 0 + 1
polg.indx <- terra::as.polygons(mask.indx)

# Extract by mask the climate 
clma <- terra::crop(clma, polg.indx) %>% terra::mask(., polg.indx)
indx <- terra::crop(indx, polg.indx) %>% terra::mask(., polg.indx)

# To create the stack
stck <- c(clma, indx)
terra::writeRaster(x = stck, filename = '../data/tif/results/climate_indx-raw_ssp370_col-ecu-per_v1.tif', overwrite = TRUE)

names(stck)

rast('../data/tif/results/indx_bsl-ftr_chirts-ssp370.tif')

