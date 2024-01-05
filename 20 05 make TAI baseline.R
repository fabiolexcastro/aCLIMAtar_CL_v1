

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
rdrs <- function(x){
  grep(x, dirs, value = T) %>% 
    dir_ls() %>% 
    grep('monthly', ., value = T) %>% 
    dir_ls() %>% 
    rast()
}
ctai <- function(y){
  
  # To start the function
  cat('... To process: ', y, '\n')
  
  # Filtering the raster for each year
  ppt <- prec[[grep(y, names(prec))]]
  tmn <- tmin[[grep(y, names(tmin))]]
  tmx <- tmax[[grep(y, names(tmax))]]
  
  rst <- map(.x = 1:12, .f = function(i){
    
    # Starting 
    cat('To process: ', month.abb[i], '\n')
    
    # Select the months
    pp <- ppt[[i]]
    tn <- tmn[[i]]
    tx <- tmx[[i]]
    
    # To calculate average temperature
    ta <- (tn + tx) / 2
    
    # To calculate the range
    rn <- abs(tx - tn)
    rn <- mean(rn, na.rm = T)
    
    # Return the object
    cat('Done!\n')
    return(list(pp, ta, rn))
    
  })
  
  # To extract the object
  ppt <- map(rst, 1)
  tav <- map(rst, 2)
  rng <- map(rst, 3)
  srd <- srad
  
  # To create the stack 
  ppt <- reduce(ppt, c)
  tav <- reduce(tav, c)
  rng <- reduce(rng, c)

  # To change the names
  names(ppt) <- c(glue('PREC_0{1:9}'), glue('PREC_{10:12}'))
  names(tav) <- c(glue('TMEAN_{1:9}'), glue('TMEAN_{10:12}'))
  names(srd) <- c(glue('SRAD_0{1:9}'), glue('SRAD_{10:12}'))
  names(rng) <- c(glue('TRNG_0{1:9}'), glue('TRNG_{10:12}'))
  
  # Assign precpitation names in enviren environment
  envirem::assignNames(solrad = 'SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
  
  # Rast to raster
  ppt <- raster::stack(ppt)
  tav <- raster::stack(tav)
  rng <- raster::stack(rng)
  srd <- raster::stack(srd)
  
  # Resampling
  srd <- raster::resample(srd, tav, method = 'bilinear')
  
  # To calculate PET
  pet <- envirem::monthlyPET(tav, srd, rng)
  names(pet)
  
  # To calculate TAI
  tai <- envirem::aridityIndexThornthwaite(ppt, pet)
  tai <- rast(tai)
  names(tai) <- glue('tai_{y}')
  
  # To return
  rm(pet, tav, rng); gc(reset = T)
  return(tai)
    
}

# Load data ---------------------------------------------------------------

# Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('ECU', 'PER', 'COL')
zone <- wrld[wrld$sov_a3 %in% isos,]
zone <- vect(zone)

# Climate data ----------------

# Solar radiation 
srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad', regexp = 'solrad_')
srad <- mixedsort(srad)
srad <- rast(srad)
srad <- terra::crop(srad, zone)
srad <- terra::mask(srad, zone)
srad <- terra::crop(srad, c(-82, -66.9, -18.4, 12.5))

# Mains vars
dirs <- dir_ls('../data/tif/climate/baseline', type = 'directory')
prec <- rdrs('prec')
tmin <- rdrs('tmin')
tmax <- rdrs('tmax')

grep('tmax', dirs, value = T) %>% dir_ls() %>% grep('monthly', ., value = T) %>% dir_ls() %>% as.character() %>% basename() %>% str_sub(., 6, 9) %>% table()

# To change the names
dtes <- seq(as.Date('1983-01-01', format = '%Y-%m-%d'), as.Date('2015-12-31', format = '%Y-%m-%d'), by = 'month')
names(tmin) <- glue('tmin_{dtes}')
names(tmax) <- glue('tmax_{dtes}')

# To calculate TAI --------------------------------------------------------

# Time series
tais <- map(1983:2015, ctai)
tais <- reduce(tais, c)
terra::writeRater(tais, '../data/tif/index/baseline/tai/tai_bsl_tsr_chirts.tif', overwrite = TRUE)

# Average
tais.avrg <- app(tais, mean)
names(tais.avrg) <- 'tai'
terra::writeRaster(x = tais.avrg, '../data/tif/index/baseline/tai/tai_bsl_avg_chirts.tif', overwrite = TRUE)

