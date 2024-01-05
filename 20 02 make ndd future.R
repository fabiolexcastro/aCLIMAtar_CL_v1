

## Make index: number of dry days (NDD) - Future dataset
## Source: https://github.com/AdaptationAtlas/hazards/wiki/Hazards-definitions
## Author: Fabio Castro - Llanos 
## Alliance Bioversity - CIAT

### Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rgeos, crayon, gtools, stringr, glue, future, furrr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

### Load data ---------------------------------------------------------------

# Spatial data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', ' ECU', 'PER')

# Raster data
dirs <- dir_ls('../data/tif/climate/future/prec/ssp370')
mdls <- basename(dirs)

### Function to use ---------------------------------------------------------
calc.ndd <- function(mdl){
  
  # mdl <- mdls[1]
  # To filtering the files
  cat('To process: ', mdl, '\n')
  fles <- dir_ls(dir_ls(as.character(grep(mdl, dirs, value = T)), regexp = 'daily'))
  year <- unique(parse_number(basename(fles)))
  fles <- as.character(grep(paste0(paste0('_', 2026:2055, '-'), collapse = '|'), fles, value = T))
  dtes <- basename(fles) %>% str_split(., pattern = '_') %>% map_chr(2) %>% gsub('.tif$', '', .)
  
  # To calculate the number of dry days for each month in each year
  nddr <- map(.x = 1:length(fles), .f = function(i){
    cat('NDD: processing: ', i, '\n')
    fle <- fles[i]
    yrr <- parse_number(basename(fle))
    mnt <- str_sub(basename(fle), 11, 12)
    rst <- rast(fle)
    ndd <- terra::app(x = rst, fun = function(x){ ndd <- sum(x < 1, na.rm = T); return(ndd)})
    return(ndd)
  })
  
  nddr <- reduce(nddr, c)
  names(nddr) <- glue('ndd_{dtes}')
  
  # To calculate the average for each month
  ndda <- map(.x = 1:12, .f = function(m){
    m <- ifelse(m < 10, paste0('0', m), as.character(m))
    r <- nddr[[grep(m, names(nddr))]]
    r <- mean(r)
    names(r) <- glue('ndd_{m}')
    return(r)
  }) 
  
  # To reduce and calculate the final raster
  ndda <- reduce(ndda, c)
  ndda <- sum(ndda)
  ndda <- floor(ndda)
  return(ndda)
  
}

### To apply the function ---------------------------------------------------
rstr <- map(.x = 1:length(mdls), .f = function(x){
  mdl <- mdls[x]
  ndd <- calc.ndd(mdl = mdl)
  return(ndd)
})
rstr <- reduce(rstr, c)
names(rstr) <- glue('ndd_{mdls}')
dir_create('../data/tif/index/future/ndd')
terra::writeRaster(x = rstr, filename = '../data/tif/index/future/ndd/ndd_ftr-ssp370-mdls_raw_all.tif', overwrite = TRUE)

rast('../data/tif/index/future/ndd/ndd_ftr-ssp370-mdls_raw_all.tif')
