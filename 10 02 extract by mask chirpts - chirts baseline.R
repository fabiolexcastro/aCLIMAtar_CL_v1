
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue, stringr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Path raster data
chrp <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps'
chrt <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global'

# Administrative data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Function to extract by mask ---------------------------------------------

# Extract by mask ---------------------------------------------------------

# CHIRPS ------------------------------------------------------------------
fles.prec <- dir_ls(chrp, regexp = '.tif') %>% as.character() 
year.prec <- basename(fles.prec) %>% str_split(pattern = '\\.') %>% map_chr(3) %>% as.numeric() %>% unique()
map(.x = year.prec, .f = function(y){
  cat('To extract by mask: ', y)
  rstr <- rast(grep(y, fles.prec, value = T))
  fnal <- map(.x = 1:12, .f = function(m){
    cat('To extract by mask: ', m, '\n')
    m <- ifelse(m < 10, paste0('0', m), as.character(m))
    r <- rstr[[grep(paste0(y, '.', m, '.'), names(rstr))]]
    r <- terra::crop(r, vect(zone))
    r <- terra::mask(r, vect(zone))
    r <- terra::crop(r, c(-82, -66.9, -18.4, 12.5))
    r <- ifel(r < 0, 0, r)
    a <- sum(r)
    return(list(r, a))
  })
  trra <- map(fnal, 1)
  mnth <- map(fnal, 2)
  
  trra <- reduce(trra, c)
  mnth <- reduce(mnth, c)
  names(mnth) <- glue('prec_{y}-{1:12}')
  
  out <- glue('../data/tif/climate/baseline/prec')
  terra::writeRaster(x = trra, filename = glue('{out}/daily/prec_{y}.tif'), overwrite = TRUE)
  terra::writeRaster(x = mnth, filename = glue('{out}/monthly/prec_{y}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
})


# Total summarisse --------------------------------------------------------
rstr <- dir_ls('../data/tif/climate/baseline/prec', regexp = 'monthly') %>% 
  dir_ls() %>% 
  as.character() %>% 
  rast()
rstr <- map(.x = 1:12, .f = function(i){mean(rstr[[grep(paste0('-', i), names(rstr))]])})
rstr <- reduce(rstr, c) 
names(rstr) <- glue('prec_{1:12}')

# CHIRTS ------------------------------------------------------------------
fles.tasm <- dir_ls(chrt, type = 'directory') %>% as.character()
vars <- c('Rh', 'Tmax', 'Tmin')

extr.chrt <- function(var){
  
  var <- 'Tmax'
  
  fls <- grep(var, fles.tasm, value = TRUE) %>% 
    dir_ls() %>% 
    as.character()
  
  yrs <- 1983:2015
  
  map(.x = 1:length(yrs), .f = function(y){
    
    # year <- 1983 # Run after
    year <- yrs[y]
    fles <- grep(paste0('.', year, '.'), fls, value = T)
    
    map(.x = 1:12, .f = function(m){
      
      cat('To extract by mask: ', year, ' ', m, '\n')
      
      # To read and extract by mask
      m <- ifelse(m < 10, paste0('0', m), as.character(m))
      r <- grep(paste0(year, '.', m, '.'), fles, value = T)
      r <- rast(r)
      r <- terra::crop(r, vect(zone))
      r <- terra::mask(r, vect(zone))
      r <- ifel(r < -9998, NA, r)
      r <- terra::crop(r, c(-82, -66.9, -18.4, 12.5))
      a <- mean(r)
      
      # To write
      o <- glue('../data/tif/climate/baseline/{tolower(var)}')
      dir_create(glue('{o}/monthly'))
      dir_create(glue('{o}/daily'))
      terra::writeRaster(x = r, filename = glue('{o}/daily/{tolower(var)}_{year}-{m}.tif'), overwrite = TRUE)
      terra::writeRaster(x = a, filename = glue('{o}/monthly/{tolower(var)}_{year}-{m}.tif'), overwrite = TRUE)
      cat('Done!\n')
      rm(r, a, o)
      gc(reset = TRUE)
      cat('Done!\n')
      
    })
    
  })
  
  
}


# RH ----------------------------------------------------------------------
fles <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/ERA5/2m_relative_humidity', regexp = '.nc$')
fles <- mixedsort(fles)
fles <- as.character(fles)
nmes <- basename(fles) %>% str_split(string = ., pattern = '_') %>% map_chr(4)
nmes %>% as.numeric() %>% mixedsort()

# Sort by each year
year <- 1981:2022
fles <- fles[order(nmes)]

map(.x = 1:length(year), .f = function(i){
  
  # To read as a raster file
  # file <- fles[1]
  cat('To process: ', year[i], '\n')
  yr  <- year[i]
  flss <- grep(paste0('_', yr), fles, value = T)
  
  rstr <- map(.x = 1:length(fls), .f = function(j){
  
    cat('To process: ', j, '\n')
    file <- flss[j]
    rstr <- rast(file)
    
    # To extract by mask
    rstr <- terra::crop(rstr, zone)
    rstr <- terra::mask(rstr, zone)
    rstr <- terra::crop(rstr, c(-82, -66.9, -18.4, 12.5))
    name <- basename(file)
    dtes <- str_split(name, '_') %>% map_chr(4)
    
    # Get the dates
    yeal <- str_sub(dtes, 1, 4)
    mnth <- str_sub(dtes, 5, 6)
    days <- str_sub(dtes, 7, 8)
    
    # To return the raster
    names(rstr) <- glue('rhum_{yeal}-{mnth}-{days}')
    return(rstr)
    
  })
  
  rstr <- reduce(rstr, c)
  plot(rstr)
  dout <- glue('../data/tif/climate/baseline/rhum/daily')
  terra::writeRaster(x = rstr, filename = glue('{dout}/rhum_{yr}.tif'), overwrite = TRUE)
  
  # Now daily to monthly
  avrg <- map(.x = 1:12, .f = function(m){
    
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    rst <- rstr[[grep(paste0('-', mnt, '-'), names(rstr))]]
    avg <- mean(rst)
    names(avg) <- glue('rhum_{yr}-{mnt}')
    return(avg)
    
  })
  
  avrg <- reduce(avrg, c)
  dout <- glue('../data/tif/climate/baseline/rhum/monthly')
  terra::writeRaster(x = avrg, filename = glue('{dout}/rhum_{yr}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
  rm(rstr, avrg); gc(reset = TRUE)
  
})
