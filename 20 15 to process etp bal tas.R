

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue, gtools, stringr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions ---------------------------------------------------------------
calc.etps <- function(prd){
  
  # Start the analysis
  cat('To process: ', prd, '\n')
  stk <- stck[[grep(prd, names(stck), value = F)]]
  ppt <- stk[[grep('prec', names(stk), value = F)]]
  tmn <- stk[[grep('tmin', names(stk), value = F)]]
  tmx <- stk[[grep('tmax', names(stk), value = F)]]
  
  # To calculate average temperature
  tav <- (tmx + tmn) / 2
  
  # To calculate the ETP
  etp <- 0.0013 * 0.408 * srad * (tav + 17) * (tmx - tmn - 0.0123 * ppt) ^ 0.76
  etp <- etp * c(31,29,31,30,31,30,31,31,30,31,30,31)
  for(i in 1:12){etp[[i]][which.lyr(is.na(etp[[i]]))] <- 0}
  etp <- terra::crop(etp, zone) %>% terra::mask(., zone)
  names(etp) <- glue('etps_{prd}_{c(paste0(0, 1:9), 10:12)}')
  
  # Return the raster stack
  cat('Done!\n')
  return(etp)
  
}
calc.baln <- function(prd){
  
  cat('To process: ', prd, '\n')
  # Evapotranspiration
  etp <- stck[[grep('etps', names(stck), value = F)]]
  etp <- etp[[grep(prd, names(etp), value = F)]]
  
  # Precipitation
  ppt <- stck[[grep('prec', names(stck), value = F)]]
  ppt <- ppt[[grep(prd, names(ppt), value = F)]]
  
  # To calculate the balance
  bal <- ppt - etp
  
  # To make the reclassify
  mtx <- matrix(c(-Inf, -10, -1, -10, 10, 0, 0, Inf, 1), ncol = 3, nrow = 3, byrow = T)
  rcl <- terra::classify(bal, mtx)
  names(rcl) <- glue('baln_{prd}_{c(paste0(0, 1:9), 10:12)}')
  
  # Finish
  cat('Done!')
  return(rcl)
  
}
calc.tavg <- function(prd){
  
  cat('To process: ', prd, '\n')
  rst <- stck[[grep(prd, names(stck))]]
  tmn <- rst[[grep('tmin', names(rst))]]
  tmx <- rst[[grep('tmax', names(rst))]]
  tav <- (tmn + tmx) / 2 
  names(tav) <- glue('tavg_{prd}_{1:12}')
  return(tav)
  
}
calc.dfrn <- function(varb){
  
  # To start the analysis - filtering the raster
  cat('To process: ', varb, '\n')
  rstr <- stck[[grep(varb, names(stck))]]
  ftre <- rstr[[grep('ftr', names(rstr), value = FALSE)]]
  bsln <- rstr[[grep('bsl', names(rstr), value = FALSE)]] 
  
  # To calculate the difference
  dfrn <- ftre - bsln
  names(dfrn) <- glue('{varb}_ftr-dfr_{1:12}')
  
  # Finish
  cat('Done!')
  return(dfrn)
  
}
rst2tbl <- function(stk){
  
  cat('To process!\t')
  tbl <- terra::as.data.frame(stk, xy = T) %>% 
    as_tibble() %>% 
    drop_na()
  
  cat('Done!\n')
  return(tbl)
  
}

# To read climate data ----------------------------------------------------

# Raster data
stck <- terra::rast('../data/tif/results/climate_indx-raw_ssp370_col-ecu-per_v1.tif')

# Vector data
zone <- vect('../data/gpk/zone.gpkg')

# Solar radiation
srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad')
srad <- mixedsort(srad)
srad <- grep('et_solrad', srad, value = T)
srad <- as.character(srad)
srad <- rast(srad)

# To extract by mask solar radiation
srad <- terra::crop(srad, zone)
srad <- terra::mask(srad, zone)
srad <- terra::crop(srad, ext(stck))
srad <- terra::resample(srad, stck)

# To calculate ETP --------------------------------------------------------
etps.bsln <- calc.etps(prd = 'bsl')
etps.frcs <- calc.etps(prd = 'frc')
etps.ftre <- calc.etps(prd = 'ftr-avg')

# Add the ETP to the stack 
stck <- c(stck, etps.bsln, etps.frcs, etps.ftre)

# To calculate baln -------------------------------------------------------
baln.bsln <- calc.baln(prd = 'bsl')
baln.frcs <- calc.baln(prd = 'frc')
baln.ftre <- calc.baln(prd = 'ftr-avg')

# Add the baln to the stack 
stck <- c(stck, baln.bsln, baln.frcs, baln.ftre)

# To calculate the average temperature ------------------------------------

# To calculate average temperature 
tavg.bsln <- calc.tavg(prd = 'bsl')
tavg.frcs <- calc.tavg(prd = 'frc')
tavg.ftre <- calc.tavg(prd = 'ftr')

# Add the tavg to the stack
stck <- c(stck, tavg.bsln, tavg.frcs, tavg.ftre)

# Temperature difference --------------------------------------------------

# To calculate
tmin.dfrn <- calc.dfrn(varb = 'tmin')
tavg.dfrn <- calc.dfrn(varb = 'tavg')
tmax.dfrn <- calc.dfrn(varb = 'tmax')

# Get the threshold
sd.tmax <- as.numeric(global(app(tmax.dfrn, 'sd'), 'mean', na.rm = T))
sd.tmin <- as.numeric(global(app(tmin.dfrn, 'sd'), 'mean', na.rm = T))
sd.tavg <- as.numeric(global(app(tavg.dfrn, 'sd'), 'mean', na.rm = T))

tasm.binr <- function(stk, thr){
  rsl <- map(.x = 1:12, .f = function(i){
    rcl <- stk[[i]]
    rcl <- terra::classify(rcl, matrix(c(-Inf, -1*thr, -1, -1*thr, sd.tmax, 0, thr, Inf, 1), byrow = T, ncol = 3))
    return(rcl)
  }) %>% 
    reduce(., c)
  return(rsl)
}

# Add the difference to the stack 
tmin.dfrn <- tasm.binr(stk = tmin.dfrn, thr = sd.tmin)
tavg.dfrn <- tasm.binr(stk = tavg.dfrn, thr = sd.tavg)
tmax.dfrn <- tasm.binr(stk = tmax.dfrn, thr = sd.tmax)

names(stck)
stck <- c(stck, tmin.dfrn, tavg.dfrn, tmax.dfrn)

# To write the final raster -----------------------------------------------
terra::writeRaster(x = stck, filename = '../data/tif/results/climate_indx-raw_ssp370_col-ecu-per_v2.tif', overwrite = TRUE)

# To read the index -------------------------------------------------------
indx.col <- terra::rast('../data/tif/results/indx_bsl-ftr_chirts-ssp370_col.tif')
indx.ecu <- terra::rast('../data/tif/results/indx_bsl-ftr_chirts-ssp370_ecu.tif')
indx.per <- terra::rast('../data/tif/results/indx_bsl-ftr_chirts-ssp370_per.tif')

indx.nmes <- c('ndd_bsl', 'hsh_bsl', 'tai_bsl', 'ntx35_bsl', 'ndwl0_bsl', 'ndd_ftr', 'hss_ftr', 'tai_ftr', 'ntx35_ftr', 'ndwl0_ftr')

# Climate stack
clma <- stck[[-grep(paste0(indx.nmes, collapse = '|'), names(stck))]]

# Index by each country

## Colombia 
clma.col <- terra::crop(clma, zone[zone$sov_a3 == 'COL',]) %>% terra::mask(., zone[zone$sov_a3 == 'COL',])
stck.col <- c(clma.col, indx.col)

## Ecuador
clma.ecu <- terra::crop(clma, zone[zone$sov_a3 == 'ECU',]) %>% terra::mask(., zone[zone$sov_a3 == 'ECU',])
stck.ecu <- c(clma.ecu, indx.ecu)

## Perú 
clma.per <- terra::crop(clma, zone[zone$sov_a3 == 'PER',]) %>% terra::mask(., zone[zone$sov_a3 == 'PER',])
stck.per <- c(clma.per, indx.per)

# Raster to table ---------------------------------------------------------
tble.col <- rst2tbl(stk = stck.col)
tble.ecu <- rst2tbl(stk = stck.ecu)
tble.per <- rst2tbl(stk = stck.per)

# To write the final table ------------------------------------------------
write.csv(tble.col, '../data/tbl/results/tble_climate-index_col.csv', row.names = FALSE)
write.csv(tble.ecu, '../data/tbl/results/tble_climate_index_ecu.csv', row.names = FALSE)
write.csv(tble.per, '../data/tbl/results/tble_climate_index_per.csv', row.names = FALSE)












