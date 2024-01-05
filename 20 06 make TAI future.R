

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
rdrs <- function(x, m){
  f <- grep(x, dirs, value = T) %>% 
    dir_ls() %>% 
    grep('monthly', ., value = T) %>% 
    dir_ls() %>% 
    grep(m, ., value = T) %>% 
    as.character() 
  r <- rast(f)
  names(r) <- gsub('.tif', '', basename(f))
  return(r)
}
ctai <- function(y, m){
  
  # To start the function
  cat('... To process: ', y, m, '\n')
  
  # y <- 2026
  # m <- 'INM-CM5-0'
  
  # Grepping the files
  ppt <- rdrs(x = 'prec', m = m)
  tmn <- rdrs(x = 'tmin', m = m)
  tmx <- rdrs(x = 'tmax', m = m)
  
  # Filtering the raster for each year
  ppt <- ppt[[grep(y, names(ppt))]]
  tmn <- tmn[[grep(y, names(tmn))]]
  tmx <- tmx[[grep(y, names(tmx))]]
  
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
  names(tai) <- glue('tai_{y}_{m}')
  
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
dirs <- dir_ls('../data/tif/climate/future', type = 'directory') %>% map(., dir_ls) %>% unlist() %>% as.character() %>% map(., dir_ls) %>% unlist() %>% as.character()
mdls <- unique(basename(dirs))

# To calculate tai --------------------------------------------------------

# All the models
tais <- map(.x = 1:length(mdls), .f = function(j){
  cat('To process: ', mdls[j], '\t')
  rst <- map2(.x = 2026:2055, .y = rep(mdls[j], length(2026:2055)), .f = ctai)
  return(rst)
})
tais <- unlist(tais)
tais <- reduce(tais, c)
mdls

map(.x = 1:5, .f = function(i){
  cat('To processing: ', mdls[i], '\n')
  rst <- tais[[grep(mdls[i], names(tais))]]
  terra::writeRaster(x = rst, filename = glue('../data/tif/index/future/tai/tai_{mdls[i]}_monthly.tif'), overwrite = TRUE)
  cat('Done!\n')
})

# To calculate the average
tais
years <- names(tais) %>% 
  str_split(., pattern = '_') %>% 
  map_chr(., pluck, 2) %>% 
  unique() 

tais.avrg <- map(.x = years, .f = function(y){
  
  cat('Processing: ', y, '\n')
  r <- tais[[grep(y, names(tais))]]
  r <- app(r, 'mean')
  names(r) <- glue('tai_{y}')
  return(r)  
  
}) %>% 
  reduce(., c)


tais.avrg <- terra::app(tais.avrg, 'mean')
terra::writeRaster(x = tais.avrg, filename = '../data/tif/index/future/tai/tai_avrg.tif', overwrite = TRUE)



