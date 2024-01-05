

## Make index: number of dry days (NDD)
## Source: https://github.com/AdaptationAtlas/hazards/wiki/Hazards-definitions
## Author: Fabio Castro - Llanos 
## Alliance Bioversity - CIAT

### Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rgeos, crayon, gtools, stringr, glue, future, furrr)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

### Load data ---------------------------------------------------------------

# Raster data
path <- '../data/tif/climate/baseline/prec/daily'
fles <- dir_ls(path, regexp = '.tif$')
fles <- as.character(fles)
year <- parse_number(basename(fles))

# Spatial data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', ' ECU', 'PER')

### Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.ndd <- function(yr){
  
  # To read as a raster file
  cat('To proces::', green(yr), '\n')
  rst <- rast(grep(yr, fles, value = T))
  
  # To calculate the number of dry days by each month
  trr <- map(.x = 1:12, .f = function(m){
    cat('To process: ', m, '\n')
    m <- mnth(m)
    r <- rst[[grep(paste0(yr, '.', m,  '.'), names(rst), value = FALSE)]]
    f <- terra::app(x = r, fun = function(x){ ndd <- sum(x < 1, na.rm = T); return(ndd)})
    return(f)  
  }) 
  
  # To calculate the sum and return the raster file
  trr <- reduce(trr, c)
  trr <- sum(trr)
  names(trr) <- glue('ndd_{yr}')
  return(trr)
  
}

### To apply the function  --------------------------------------------------
nddr <- map(.x = 1990:2022, .f = calc.ndd)
nddr <- reduce(nddr, c)
nddr

# Output directory and write the raster
dout <- '../data/tif/index/baseline/ndd'
dir_create(dout)
terra::writeRaster(x = nddr, filename = glue('{dout}/ndd_bsl_raw_all.tif'), overwrite = TRUE)

# To calculate the average ------------------------------------------------
ndda <- terra::app(nddr, mean)
ndda <- terra::app(ndda, function(x) round(x, digits = 0))
names(ndda) <- glue('ndd_bsl')
terra::writeRaster(x = ndda, filename = glue('{dout}/ndd_bsl_avg_all.tif'), overwrite = TRUE)


var <- 'tmax'
deparse(quote(var))
