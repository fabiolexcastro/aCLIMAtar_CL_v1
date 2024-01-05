

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
cntx <- function(d, y, m, thr = 35){
  
  # d <- tmax.dirs[1]
  # y <- 2026
  # m <- 1
  
  # To start the analysis
  cat('To processing: ', basename(d), y, m, '\t')
  
  # Filtering for the year and month
  mnt <- mnth(m)
  drs <- as.character(dir_ls(tmax.dirs, regexp = 'daily'))
  dir <- grep(d, drs, value = T)
  fls <- dir_ls(dir, regexp = '.tif$')
  fls <- as.character(fls)
  fle <- grep(paste0('_', y, '-'), fls, value = T) 
  fle <- grep(paste0('-', mnt, '.tif'), fle, value = T)
  
  # Read as a raster file
  rst <- terra::rast(fle)
  
  # To calculate the index
  rsl <- 
    terra::app(x = rst, 
               fun = function(x){
                 ntxval = sum(x >= thr, na.rm = T)
                 return(ntxval)
               }) %>% 
    terra::crop(
      ., 
      zone
    ) %>% 
    terra::mask(
      ., 
      zone
    )
  
  names(rsl) <- glue('ntxr35_{y}-{m}')
  
  # To return the file
  cat('Done\n')
  return(rsl)
  
}
cavr <- function(stck){
  
  cat('To start the analysis\n')
  
  rst <- map(.x = 1:12, .f = function(m){
    
    cat(month.abb[m], '\t')
    stk <- stck[[grep(paste0('-', m, '$'), names(stck))]]
    avg <- terra::app(stk, 'mean')
    names(avg) <- glue('ntx35_{m}')
    return(avg)
    
  }) %>% 
    reduce(., c)
  
  return(rst)
  
}

# Load data ---------------------------------------------------------------

# Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('ECU', 'PER', 'COL')
zone <- wrld[wrld$sov_a3 %in% isos,]
zone <- vect(zone)

# Climate data ----------------
tmax.dirs <- dir_ls('../data/tif/climate/future/tmax/ssp370') %>% as.character() 

# To calculate ntx --------------------------------------------------------
years <- 2026:2055
month <- 1:12

# To apply for each model
map(.x = 1:5, .f = function(d){
  
  mdl <- tmax.dirs[d] %>% 
    basename()
  
  cat(mdl, '\t')
  
  rstr <- map(.x = 1:length(years), .f = function(i){
    
    rst <- map(.x = 1:12, .f = function(j){
      r <- cntx(d = tmax.dirs[d], y = years[i], m = month[j])
    }) %>% 
      reduce(., c)
    
  }) %>% 
    reduce(., c)
  
  terra::writeRaster(x = rstr, filename = glue('../data/tif/index/future/ntx35/ntx35_{mdl}_monthly.tif'), overwrite = TRUE)
  
})

# To read the results -----------------------------------------------------
rstr <- dir_ls('../data/tif/index/future/ntx35', regexp = '.tif$') %>% as.character() %>% map(.x = ., .f = rast)
rstr.avrg <- map(rstr, cavr)
rstr.avrg <- reduce(rstr.avrg, c)
rstr.avrg <- map(.x = 1:12, .f = function(m){
  
  rst <- rstr.avrg[[grep(paste0('_', m, '$'), names(rstr.avrg))]]
  avg <- app(avg, 'mean')
  names(avg) <- glue('ntx35_{m}')
  return(avg)
  
})
rstr.avrg <- reduce(rstr.avrg, c)
rstr.avrg <- app(rstr.avrg, 'mean')
names(rstr.avrg) <- 'ntx35_ftr'

# To write the raster 
terra::writeRaster(x = rstr.avrg, filename = '../data/tif/index/future/ntx35/ntx35_avrg_monthly.tif', overwrite = TRUE)
