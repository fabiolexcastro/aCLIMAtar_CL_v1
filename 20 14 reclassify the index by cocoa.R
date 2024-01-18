

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -2)

# Load data ---------------------------------------------------------------

# Raster data
stck <- terra::rast('../data/tif/results/climate_indx-raw_ssp370_col-ecu-per_v2.tif')

# Vector data
zone <- vect('../data/gpk/zone.gpkg')

# Tabular data
pnts.col <- read_csv('../data/tbl/presences/occ_col.csv')
pnts.ecu <- vect('../data/shp/points/points_ecuador.shp')
pnts.per <- read_csv('../data/tbl/presences/occ_per.csv')

plot(zone)
points(pnts.col$lon, pnts.col$lat, pch = 16, col = 'blue')
plot(pnts.ecu, col = 'green', pch = 16, add = T)
points(pnts.per$x, pnts.per$y, pch = 16, col = 'red')

# Extract the index from the rasters --------------------------------------
nmes.indx <- names(stck)[133:nlyr(stck)]
indx <- stck[[nmes.indx]]
mask <- stck[[1]] * 0 + 1

# To cleaning the points --------------------------------------------------

# Colombia
cll.col <- terra::extract(mask, pnts.col[,c(1, 2)], cell = TRUE)$cell
pnts.col <- pnts.col[!duplicated(cll.col),]
pnts.col <- dplyr::select(pnts.col, lon, lat)
pnts.col <- mutate(pnts.col, country = 'COL')

# Ecuador 
cll.ecu <- terra::extract(mask, pnts.ecu, cell = TRUE)$cell
pnts.ecu <- pnts.ecu[!duplicated(cll.ecu),]
pnts.ecu <- crds(pnts.ecu)
pnts.ecu <- as_tibble(pnts.ecu)
pnts.ecu <- setNames(pnts.ecu, c('lon', 'lat'))
pnts.ecu <- mutate(pnts.ecu, country = 'ECU')

# Perú
cll.per <- terra::extract(mask, pnts.per[,c(1, 2)], cell = TRUE)$cell
pnts.per <- pnts.per[!duplicated(cll.per),]
pnts.per <- dplyr::select(pnts.per, lon = x, lat = y)
pnts.per <- mutate(pnts.per, country = 'PER')

## Join all the tables into only one
pnts <- rbind(pnts.col, pnts.ecu, pnts.per)

# To make the reclassify --------------------------------------------------
isos <- c('COL', 'ECU', 'PER')

## Function -----------------------
get.rcl <- function(rst, iso, max){
  
  cat('To process: ', iso, '\n')
  # To filtering
  pnt <- filter(pnts, country == iso)
  zne <- zone[zone$sov_a3 == iso,]
  
  # To extract by mask 
  rst <- terra::crop(rst, zne)
  rst <- terra::mask(rst, zne)
  
  # To extract the values for the stack
  vls <- terra::extract(x = rst, y = pnt[,c('lon', 'lat')])
  vls <- cbind(pnt, vls)
  vls <- as_tibble(vls) 
  vls <- dplyr::select(vls, -ID)
  colnames(vls)[4:5] <- c('bsl', 'ftr')
  vls <- drop_na(vls)
  
  # To calculate the percentile
  qnt <- quantile(vls$bsl, c(0, 0.5, 0.75, 1))
  mtr <- matrix(c(0, qnt[2], 1, qnt[2], qnt[3], 2, qnt[3], max, 3), byrow = T, ncol = 3)
  
  # To reclassify
  cls <- terra::classify(rst, mtr, include.lowest = TRUE)
  
  # Finish
  cat('Done!\n')
  return(cls)
  
}

### Number of dry days ------------
ndd.raw <- indx[[c('ndd_bsl', 'ndd_ftr')]]
ndd.rcl <- map(1:length(isos), function(i){get.rcl(rst = ndd.raw, iso = isos[i], max = 366)})

### Thornwaite Aridity Index ------
tai.raw <- indx[[c('tai_bsl', 'tai_ftr')]]
tai.rcl <- map(1:length(isos), function(i){get.rcl(rst = tai.raw, iso = isos[i], max = 101)})

### NTX 35 ------------------------
ntx.raw <- indx[[c('ntx35_bsl', 'ntx35_ftr')]]
ntx.rcl <- map(1:length(isos), function(i){get.rcl(rst = ntx.raw, iso = isos[i], max = 366)})

### Human Heat Stress ------------
hsh.raw <- indx[[c('hsh_bsl', 'hss_ftr')]]
mtx <- matrix(c(-Inf, 25, 1, 25, 30, 2, 30, 35, 3, 35, Inf, 4), byrow = T, ncol = 3)
hsh.cls <- terra::classify(hsh.raw, mtx, include.lowest = T)
hsh.cls <- map(1:length(isos), function(i){
  zne <- zone[zone$sov_a3 == isos[i]]
  rst <- terra::crop(hsh.cls, zne) %>% terra::mask(., zne)
})

### NDWL0 -----------------------
ndw.raw <- indx[[c('ndwl0_bsl', 'ndwl0_ftr')]]
ndw.cls <- map(1:length(isos), function(i){get.rcl(rst = ndw.raw, iso = isos[i], max = 366)})

# To make the stack for each country --------------------------------------
ind.col <- list(ndd.rcl[[1]], tai.rcl[[1]], ntx.rcl[[1]], hsh.cls[[1]], ndw.cls[[1]])
ind.col <- reduce(ind.col, c)
ind.ecu <- list(ndd.rcl[[2]], tai.rcl[[2]], ntx.rcl[[2]], hsh.cls[[2]], ndw.cls[[2]])
ind.ecu <- reduce(ind.ecu, c)
ind.per <- list(ndd.rcl[[3]], tai.rcl[[3]], ntx.rcl[[3]], hsh.cls[[3]], ndw.cls[[3]])
ind.per <- reduce(ind.per, c)

# To write 
terra::writeRaster(x = ind.col, filename = '../data/tif/results/indx_bsl-ftr_chirts-ssp370_col.tif', overwrite = TRUE)
terra::writeRaster(x = ind.ecu, filename = '../data/tif/results/indx_bsl-ftr_chirts-ssp370_ecu.tif', overwrite = TRUE)
terra::writeRaster(x = ind.per, filename = '../data/tif/results/indx_bsl-ftr_chirts-ssp370_per.tif', overwrite = TRUE)





