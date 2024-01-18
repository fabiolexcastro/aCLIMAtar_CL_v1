
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
calc.bsln.avrg.mnth <- function(fles, varb){
  
  cat('To process: ', varb, '\n')
  fles <- grep(varb, fles, value = T)
  
  rstr <- map(.x = 1:12, .f = function(m){
    
    cat('To process: ', m, '\t')
    mnt <- ifelse(m < 10, glue('0{m}'), m)
    fls <- grep(glue('-{mnt}.tif'), fles, value = T)  
    fls <- grep(paste0(1995:2015, collapse = '|'), fls, value = T)
    rst <- terra::rast(fls)
    avg <- terra::app(rst, mean)
    names(avg) <- glue('{varb}_{mnt}')
    return(avg)
    
  }) %>% 
    reduce(., c)
  
  cat('Finish\n')
  return(rstr)
  
}
calc.ftre.avrg.mnth <- function(dirs, varb){
  
  # dirs <- dirs.ftre
  # varb <- 'tmax'
  
  cat('To process: ', varb, '\n')
  
  dirs <- grep(varb, dirs.ftre, value = T)
  fles <- paste0(dirs, '/monthly') %>% dir_ls() %>% unlist() %>% as.character()
  
  rstr <- map(.x = 1:5, .f = function(i){
    
    cat(mdls[i], '\n')
    
    dta <- map(.x = 1:12, .f = function(j){
      
      cat(month.abb[j], '\t')
      mnt <- ifelse(j < 10, paste0('0', j), j)
      fls <- grep(paste0('-', mnt, '.tif$'), fles, value = T) %>% grep(mdls[i], ., value = T) 
      rst <- terra::rast(fls)
      rst <- terra::app(rst, mean)
      names(rst) <- glue('{varb}_{mdls[m]}_{j}')
      return(rst)
      
    }) %>% 
      reduce(., c)
    
    cat('Done!\n')
    return(dta)
    
  }) %>% 
    reduce(., c)
  
  # To calculate min - avg - max
  stck <- map(.x = 1:12, .f = function(k){
    
    cat('---> Processing: ', month.abb[k], '\n')
    rst <- rstr[[grep(paste0('_', k, '$'), names(rstr))]]
    # stk <- c(min(rst), mean(rst), max(rst))
    stk <- mean(rst)
    names(stk) <- glue('{varb}_mean_{k}')
    # names(stk) <- c(glue('{varb}_min_{k}'), glue('{varb}_mean_{k}'), glue('{varb}_max_{k}'))
    return(stk)
    
  }) 
  
  
  stck <- reduce(stck, c)
  return(stck)
  
  
}
calc.ftre.prec.mnth <- function(dirs, varb){
  
  # dirs <- dirs.ftre
  # varb <- 'prec'  
  
  cat('To process: ', varb, '\n')
  dirs <- grep(varb, dirs.ftre, value = T)
  fles <- paste0(dirs, '/monthly') %>% dir_ls() %>% unlist() %>% as.character()
  
  rstr <- map(.x = 1:5, .f = function(i){
    
    cat(mdls[i], '\n')
    
    dta <- map(.x = 1:12, .f = function(j){
      
      cat(month.abb[j], '\t')
      mnt <- ifelse(j < 10, paste0('0', j), j)
      fls <- grep(paste0('-', mnt, '.tif$'), fles, value = T) %>% grep(mdls[i], ., value = T)
      rst <- terra::rast(fls)
      rst <- terra::app(rst, mean)
      names(rst) <- glue('{varb}_{mdls[m]}_{j}')
      return(rst)
      
    }) %>% 
      reduce(., c)
    
    cat('Done!\n')
    return(dta)
    
  }) %>% 
    reduce(., c)
  
  # To calculate min - avg - max
  stck <- map(.x = 1:12, .f = function(k){
    
    cat('---> Processing: ', month.abb[k], '\n')
    rst <- rstr[[grep(paste0('_', k, '$'), names(rstr))]]
    stk <- c(min(rst), mean(rst), max(rst))
    names(stk) <- c(glue('{varb}_min_{k}'), glue('{varb}_mean_{k}'), glue('{varb}_max_{k}'))
    return(stk)
    
  }) 
  
  stck <- reduce(stck, c)
  return(stck)

}


# Baseline data -----------------------------------------------------------

fles.bsln <- dir_ls('../data/tif/climate/baseline') %>% 
  paste0(., '/monthly') %>% 
  grep(paste0(c('prec', 'tmin', 'tmax'), collapse = '|'), ., value = T) %>% 
  dir_ls() %>% 
  as.character()

## Temperature -------
tmin.bsln <- calc.bsln.avrg.mnth(fles = fles.bsln, varb = 'tmin')
tmax.bsln <- calc.bsln.avrg.mnth(fles = fles.bsln, varb = 'tmax')

names(tmin.bsln) <- glue('tmin_bsl_{1:12}')
names(tmax.bsln) <- glue('tmax_bsl_{1:12}')

## Precipitation -----
prec.bsln <- grep('prec', fles.bsln, value = T) %>% grep(paste0(1995:2015, collapse = '|'), ., value = T) %>% terra::rast()
prec.bsln <- map(.x = 1:12, .f = function(m){mean(prec.bsln[[grep(paste0('-', m, '$'), names(prec.bsln), value = F)]])})
prec.bsln <- reduce(prec.bsln, c)
names(prec.bsln) <- glue('prec_bsl_{1:12}')

# Forecast data -----------------------------------------------------------
fles.frcs <- dir_ls('../data/tif/climate/forecast', regexp = '.tif$')
fles.frcs <- as.character(fles.frcs)
prec.frcs <- rast(grep('prec', fles.frcs, value = T))
tmin.frcs <- rast(grep('tmin', fles.frcs, value = T))
tmax.frcs <- rast(grep('tmax', fles.frcs, value = T))

# To change the names
names(prec.frcs) <- glue('prec_frc_{1:12}')
names(tmin.frcs) <- glue('tmin_frc_{1:12}')
names(tmax.frcs) <- glue('tmax_frc_{1:12}')

# Future data -------------------------------------------------------------
mdls <- c('ACCESS', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
dirs.ftre <- dir_ls('../data/tif/climate/future') %>% paste0(., '/ssp370') %>% map(., dir_ls) %>% unlist() %>% as.character()

## To calculate the statistics for each variable --------
# Temperature
tmin.ftre <- calc.ftre.avrg.mnth(dirs = dirs.ftre, varb = 'tmin')
tmax.ftre <- calc.ftre.avrg.mnth(dirs = dirs.ftre, varb = 'tmax')

# Precipitation
prec.ftre <- calc.ftre.prec.mnth(dirs = dirs.ftre, varb = 'prec')

# To change the names 
names(tmin.ftre) <- glue('tmin_ftr-avg_{1:12}')
names(tmax.ftre) <- glue('tmax_ftr-avg_{1:12}')

nms <- list()
for(i in 1:12){
  nms[[i]] <- c(paste0('prec_ftr-min_', i), paste0('prec_ftr-avg_', i), paste0('prec_ftr-max_', i))
}
nms <- unlist(nms)
names(prec.ftre) <- nms

# To write the results 
terra::writeRaster(x = prec.ftre, filename = '../data/tif/climate/future/prec/prec_ssp370_min-avrg-max_monthly.tif', overwrite = TRUE)
terra::writeRaster(x = tmin.ftre, filename = '../data/tif/climate/future/tmin/tmin_ssp370_avrg_monthly.tif', overwrite = TRUE)
terra::writeRaster(x = tmax.ftre, filename = '../data/tif/climate/future/tmax/tmax_ssp370_avrg_monthly.tif', overwrite = TRUE)

# Join baseline - forecast - future -----------------------------------------

prec.bsln <- terra::crop(prec.bsln, ext(prec.frcs))
prec.frcs <- terra::crop(prec.frcs, ext(prec.frcs))
prec.ftre <- terra::crop(prec.ftre, ext(prec.frcs))

tmin.bsln <- terra::crop(tmin.bsln, ext(tmin.frcs))
tmin.frcs <- terra::crop(tmin.frcs, ext(tmin.frcs))
tmin.ftre <- terra::crop(tmin.ftre, ext(tmin.frcs))

tmax.bsln <- terra::crop(tmax.bsln, ext(tmax.frcs))
tmax.frcs <- terra::crop(tmax.frcs, ext(tmax.frcs))
tmax.ftre <- terra::crop(tmax.ftre, ext(tmax.frcs))

clma <- c(prec.bsln, prec.frcs, prec.ftre, 
          tmin.bsln, tmin.frcs, tmin.ftre, 
          tmax.bsln, tmax.frcs, tmax.ftre)

terra::writeRaster(x = clma, filename = '../data/tif/climate/climate_stack_col-ecu-per.tif', overwrite = TRUE)

plot(c(prec.bsln[[1]], prec.frcs[[1]], prec.ftre[[2]]))


