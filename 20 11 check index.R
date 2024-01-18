
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
calc.avrg.mnth <- function(stk, mdl){
  
  cat('Processing...', mdl, '\t')
  prd <- 2026:2055
  stk <- stk[[grep(paste0(paste0('_', prd), collapse = '|'), names(stk))]]
  rst <- map(prd, function(y){mean(stk[[grep(paste0('_', y, '_'), names(stk), value = F)]])})
  rst <- reduce(rst, c)
  names(rst) <- glue('{mdl}_{prd}')
  avg <- app(rst, mean)
  names(avg) <- glue('{mdl}')
  cat('Done!\n')
  return(avg)
  
}
calc.avrg.year <- function(stk, mdl){
  
  cat('Processing...', mdl, '\t')
  stk <- terra::app(stk, mean)
  names(stk) <- glue('{mdl}_ftr')
  return(stk)
  
}
calc.suma.mnth <- function(stk, mdl){
  
  cat('Processing...', mdl, '\t')
  prd <- 2026:2055
  stk <- stk[[grep(paste0(paste0('_', prd), collapse = '|'), names(stk))]]
  rst <- map(prd, function(y){sum(stk[[grep(paste0('_', y, '-'), names(stk), value = F)]])})
  rst <- reduce(rst, c)
  names(rst) <- glue('{mdl}_{prd}')
  avg <- app(rst, mean)
  names(avg) <- glue('{mdl}_ftr')
  return(avg)
  
}

# Load data ---------------------------------------------------------------
path <- '../data/tif/index'
dirs <- dir_ls(path, type = 'directory')

# Baseline ----------------------------------------------------------------
fles.bsln <- grep('baseline', dirs, value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  map(., dir_ls) %>% 
  unlist() %>% 
  as.character()

## Number of Dry Days ---------------
ndd.bsl <- grep('ndd', fles.bsln, value = T) %>% 
  grep('raw_all', ., value = T) %>% 
  terra::rast(.) %>% 
  .[[grep(paste0(1995:2015, collapse = '|'), names(.), value = F)]] %>% 
  terra::app(., mean) %>% 
  terra::app(., function(x) round(x, 0)) %>% 
  setNames('ndd_bsl')

## Human Heat Stress ---------------
hsh.bsl <- grep('hsh', fles.bsln, value = T) %>% 
  terra::rast() %>% 
  .[[grep(paste0(1995:2015, collapse = '|'), names(.), value = F)]] 
hsh.bsl <- map(1:12, function(x){
    mean(hsh.bsl[[grep(paste0('_', x, '$'), names(hsh.bsl), value = F)]])
  }) %>% 
  reduce(., c) %>% 
  terra::app(., mean) %>% 
  setNames('hsh_bsl')

## TAI ---------------
tai.bsl <- grep('tai', fles.bsln, value = T) %>% 
  terra::rast() %>% 
  .[[grep(paste0(1995:2015, collapse = '|'), names(.), value = F)]] %>% 
  terra::app(., mean) %>%
  round(., 0)
names(tai.bsl) <- 'tai_bsl'

## NTX35 ---------------
ntx35.bsl <- grep('ntx35', fles.bsln, value = T) %>% 
  grep('chirts', ., value = T) %>% 
  grep('tsr', ., value = T) %>% 
  terra::rast() %>% 
  .[[grep(paste0(1995:2015, collapse = '|'), names(.), value = F)]] 
ntx35.bsl <- map(1995:2015, function(y){sum(ntx35.bsl[[grep(paste0('_', y, '-'), names(ntx35.bsl), value = F)]])}) %>% 
  reduce(., c) %>% 
  mean() %>% 
  setNames('ntx35_bsl')

## NDWL0 ---------------
ndwl0.bsl <- grep('ndwl0', fles.bsln, value = T) %>% 
  terra::rast() %>% 
  .[[grep(paste0(1995:2015, collapse = '|'), names(.), value = F)]] %>% 
  terra::app(., mean) %>% 
  setNames('ndwl0_bsl')

### To create a stack for baseline ------------------------------------------
indx.bsl <- list(ndd.bsl, hsh.bsl, tai.bsl, ntx35.bsl, ndwl0.bsl) %>% reduce(., c)
names(indx.bsl)
plot(indx.bsl)

# Future ------------------------------------------------------------------
dirs.ftre <- grep('future', dirs, value = T) %>% 
  as.character() %>% 
  dir_ls() %>% 
  dir_ls() %>% 
  as.character()

mdls <- c('ACCESS', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
idxs <- c('ndd', 'hsh', 'tai', 'ntx35', 'ndwl0')

## Number of dry days ---------------
ndd.ftr <- dirs.ftre %>% 
  grep('ndd', ., value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  terra::rast() %>% 
  terra::app(., mean) %>% 
  setNames('ndd_ftr')

## Human Heat Stress ----------------
hsh.ftr <- dirs.ftre %>% 
  grep('hss', ., value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  grep(paste0(mdls, collapse = '|'), ., value = T) %>% 
  map(., rast) %>% 
  map2(., mdls, calc.avrg.mnth) %>% 
  reduce(., c) %>% 
  mean() %>% 
  setNames('hss_ftr')

hsh.ftr %>% map2(., mdls, calc.avrg.mnth)

## TAI ------------
tai.ftr <- dirs.ftre %>% 
  grep('tai', ., value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  grep(paste0(mdls, collapse = '|'), ., value = T) %>% 
  map(., rast) %>% 
  map2(., mdls, calc.avrg.year) %>% 
  reduce(., c) %>% 
  mean() %>% 
  setNames(c('tai_ftr'))

## NTX35 --------------
ntx35.ftr <- dirs.ftre %>% 
  grep('ntx35', ., value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  grep(paste0(mdls, collapse = '|'), ., value = T) %>% 
  map(., rast) %>% 
  map2(., mdls, calc.suma.mnth) %>% 
  reduce(., c) %>% 
  mean() %>% 
  round(., 0) %>% 
  setNames(c('ntx35_ftr'))

## NDWL0 --------------
ndwl0.ftr <- dirs.ftre %>% 
  grep('ndwl0', ., value = T) %>% 
  dir_ls() %>% 
  as.character() %>% 
  grep(paste0(mdls, collapse = '|'), ., value = T) %>% 
  map(., rast) %>% 
  map2(., mdls, calc.avrg.year) %>% 
  reduce(., c) %>% 
  mean() %>% 
  setNames('ndwl0_ftr')
  
### To create a stack for future ------------
indx.ftr <- list(ndd.ftr, hsh.ftr, tai.ftr, ntx35.ftr, ndwl0.ftr) %>% reduce(., c)
plot(indx.ftr)

# Join all the rasters into only one stack --------------------------------

indx <- c(indx.bsl, indx.ftr)
indx
dir_create('../data/tif/results')
terra::writeRaster(x = indx, filename = '../data/tif/results/indx_bsl-ftr_chirts-ssp370.tif', overwrite = TRUE)




