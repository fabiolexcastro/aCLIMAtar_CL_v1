
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue, stringr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Path raster data
chrp <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirps_cmip6_America'
chrt <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_cmip6_America'

# Administrative data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Function ----------------------------------------------------------------
extr.mask.prec <- function(path, sspe, prdo, name){
  
  # path <- chrp; sspe <- 'ssp370'; prdo <- '2021_2040'; name <- 'prec'
  cat('To process: ', sspe, ' ', prdo, ' ', name, '\n')
  dirs <- path %>% dir_ls(., type = 'directory') %>% grep(sspe, ., value = T) %>% as.character() %>% grep(prdo, ., value = T)
  year <- str_split(prdo, pattern = '_') %>% unlist() %>% as.numeric()
  year <- year[1]:year[2]
  varb <- deparse(substitute(chrp))

  map(.x = dirs, .f = function(d){
    
    cat('To process: ', basename(d))
    fls <- d %>% dir_ls() %>% as.character()
    mdl <- str_split(basename(d), '_') %>% map_chr(2)
    map(.x = year, .f = function(y){ 
      
      cat('Year: ', y, '\n')
      fl <- fls[grep(paste0('0.', y, '.'), basename(fls), value = F)]
      
      map(.x = 1:12, .f = function(m){
          
        cat('To extract by mask: ', m, '\n')
        
        # To read and extract by mask
        m <- ifelse(m < 10, paste0('0', m), as.character(m))
        f <- grep(paste0('.', y, '.', m, '.'), fl, value = T)
        r <- rast(f)
        r <- terra::crop(r, vect(zone))
        r <- terra::mask(r, vect(zone))
        r <- terra::crop(r, c(-82, -66.9, -18.4, 12.5))
        s <- sum(r)
        
        # To write
        o <- glue('../data/tif/climate/future/{name}/{sspe}/{mdl}')
        dir_create(glue('{o}/monthly'))
        dir_create(glue('{o}/daily'))
        terra::writeRaster(x = r, filename = glue('{o}/daily/{name}_{y}-{m}.tif'), overwrite = TRUE)
        terra::writeRaster(x = s, filename = glue('{o}/monthly/{name}_{y}-{m}.tif'), overwrite = TRUE)
        cat('Done!\n')
        rm(r, s, o)
        gc(reset = TRUE)
        cat('Done!\n')
        
      })
      
    })

  })
  
}

extr.mask.tasm <- function(path, sspe, prod, name){
  
  path <- chrt; sspe <- 'ssp370'; prdo <- '2021_2040'; name <- 'tmax'
  cat('To process: ', sspe, ' ', prdo, ' ', name, '\n')
  dirs <- path %>% dir_ls(., type = 'directory') %>% grep(sspe, ., value = T) %>% as.character() %>% grep(prdo, ., value = T)
  year <- str_split(prdo, pattern = '_') %>% unlist() %>% as.numeric()
  year <- year[1]:year[2]
  varb <- deparse(substitute(chrp))
  dirs <- grep(str_to_title(name), dirs, value = T)
  
  map(.x = dirs, .f = function(d){
    
    cat('To process: ', d, '\n')
    fls <- d %>% dir_ls() %>% map(., dir_ls) %>% unlist() %>% as.character()
    if(name == 'tmax'){
      varb <- 'Max'
    } else { 
      varb <- 'Min'
    }
    fls <- fls %>% grep(varb, ., value = T)
    mdl <- basename(d)
    mdl <- str_split(mdl, '_')
    mdl <- map_chr(mdl, 2)
    
    map(.x = year, .f = function(y){
      
      fl <- grep(paste0('/', y, '/'), fls, value = T) %>% grep('.nc$', ., value = T)
      
      map(.x = 1:12, .f = function(m){
        
        m <- ifelse(m < 10, paste0('0', m), as.character(m))
        r <- paste0('_', y, m) %>% grep(., fl, value = T) %>% rast()
        r <- terra::crop(r, zone)
        r <- terra::mask(r, zone)
        r <- terra::crop(r, c(-82, -66.9, -18.4, 12.5))
        a <- terra::app(r, 'mean')
        
        # To write
        o <- glue('../data/tif/climate/future/{name}/{sspe}/{mdl}')
        dir_create(glue('{o}/monthly'))
        dir_create(glue('{o}/daily'))
        terra::writeRaster(x = r, filename = glue('{o}/daily/{name}_{y}-{m}.tif'), overwrite = TRUE)
        terra::writeRaster(x = a, filename = glue('{o}/monthly/{name}_{y}-{m}.tif'), overwrite = TRUE)
        cat('Done!\n')
        rm(r, a, o)
        gc(reset = TRUE)
        cat('Done!\n')
        
        
      })
      
      
    })
    
  })
  
  
}

# To extract by mask ------------------------------------------------------

# CHIRPS
extr.mask.prec(path = chrp, sspe = 'ssp370', prdo = '2021_2040', name = 'prec')
extr.mask.prec(path = chrp, sspe = 'ssp370', prdo = '2041_2060', name = 'prec')

# CHIRTS
extr.mask.tasm(path = chrt, sspe = 'ssp370', prdo = '2021_2040', name = 'tmax')
extr.mask.tasm(path = chrt, sspe = 'ssp370', prdo = '2041_2060', name = 'tmax')

extr.mask.tasm(path = chrt, sspe = 'ssp370', prdo = '2021_2040', name = 'tmin')
extr.mask.tasm(path = chrt, sspe = 'ssp370', prdo = '2041_2060', name = 'tmin')



