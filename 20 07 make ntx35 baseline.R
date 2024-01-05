

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
cntx <- function(y, m, thr = 35){
  
  # To start the analysis
  cat('To processing: ', y, m, '\t')
  
  # Filtering for the year and month
  mnt <- mnth(m)
  fle <- grep(paste0('_', y, '-'), tmax, value = T) 
  fle <- grep(paste0('-', mnt, '.tif'), fle, value = T)
  
  # Read as a raster file
  rst <- terra::rast(fle)
  
  # To calculate the index
  rsl <- 
    terra::app(x = rst, 
             fun = function(x){
               ntxval = sum(x >= thr, na.rm = T)
               return(ntxval)
             })
  
  # To return the file
  cat('Done\n')
  return(rsl)
  
}

# Load data ---------------------------------------------------------------

# Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('ECU', 'PER', 'COL')
zone <- wrld[wrld$sov_a3 %in% isos,]
zone <- vect(zone)

# Climate data ----------------
tmax <- dir_ls('../data/tif/climate/baseline/tmax/daily') %>% as.character() 

# To apply the function ---------------------------------------------------

# Time series 
ntxr <- purrr::map(.x = 1983:2015, .f = function(yr){
  rs <- map(.x = 1:12, .f = function(mn){
    r <- cntx(y = yr, m = mn, thr = 35)
    names(r) <- glue('ntxr_{yr}-{mn}')
    return(r)
  })
  rs <- reduce(rs, c)
  return(rs)
})
ntxr <- reduce(ntxr, c)
dir_ls('../data/tif/index/baseline/hsh/hsh_bsl_tsr_chirts.tif')
terra::writeRaster(x = ntxr, filename = '../data/tif/index/baseline/ntx35/ntx35_bsl_tsr_chirts.tif', overwrite = TRUE)

# To calculate the average
ntxr.avrg <- ntxr[[grep(paste0(1995:2015, collapse = '|'), names(ntxr))]]
ntxr.avrg <- terra::app(ntxr.avrg, mean)
ntxr.avrg <- terra::crop(ntxr.avrg, zone)
ntxr.avrg <- terra::mask(ntxr.avrg, zone)
names(ntxr.avrg) <- 'ntxr35_bsl_avg'
terra::writeRaster(x = ntxr.avrg, filename = '../data/tif/index/baseline/ntx35/ntx35_bsl_avg_chirts.tif', overwrite = TRUE)

# -------------------------------------------------------------------------



