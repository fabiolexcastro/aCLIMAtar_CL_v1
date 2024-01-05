
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function ----------------------------------------------------------------
peest <- function(srad, tmin, tmean, tmax){
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  
  # Soil heat flux parameters
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Net radiation
  rn <- (1-albedo) * srad
  
  # Soil heat flux
  eslope <- a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
  
  # Estimate vpd
  esat_min <- 0.61120*exp((17.67*tmin)/(tmin+243.5))
  esat_max <- 0.61120*exp((17.67*tmax)/(tmax+243.5))
  vpd <- vpd_cte*(esat_max-esat_min) #kPa
  
  # Priestley-Taylor
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  pt_coef <- pt_fact*pt_const
  pt_coef <- 1 + (pt_coef-1) * vpd / vpd_ref
  
  #*10^6? To convert fluxes MJ to J
  #rlat_ht? Latent heat flux to water flux
  #100/rho_w? Kg/m^2 to cm
  et_max <- (pt_coef * rn * eslope/(eslope+psycho) * 10^6 / rlat_ht * 100/rho_w)*10 #in mm
  return(et_max)
}
eabyep_calc <<- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
  
  avail   <- min(avail, soilcp)
  
  # ERATIO
  percwt <- min(avail/soilcp*100, 100)
  percwt <- max(percwt, 1)
  eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
  
  demand  <- eratio * evap
  result  <- avail + rain - demand
  logging <- result - soilcp
  logging <- max(logging, 0)
  logging <- min(logging, soilsat)
  # runoff  <- result - logging + soilcp
  avail   <- min(soilcp, result)
  avail   <- max(avail, 0)
  # runoff  <- max(runoff, 0)
  
  out     <- list(Availability = c(AVAIL, avail),
                  # Demand       = demand,
                  # Eratio       = eratio,
                  Logging      = logging
  )
  return(out)
}
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.ndwl0 <- function(md, yr, mn){
  
  # md <- mdls[1]
  # yr <- 2028
  # mn <- 2
  
  cat('To process: ', yr, ' ', mn, '\n')
  mnt <- mnth(mn)
  
  # List the files
  drs <- dir_ls(root) %>% dir_ls() %>% as.character() %>% grep(md, ., value = T)
  fls <- map(drs, dir_ls, regexp = 'daily') %>% unlist() %>% dir_ls(., regexp = '.tif$') %>% as.character()
  
  # Period 
  yr.ft <- yr
  yr.bs <- filter(prds, future == yr) %>% pull(., 1)
  
  # Solar radiation
  srd <- grep(yr.bs, fles.srad, value = T) %>% 
    dir_ls() %>% 
    grep(paste0('-', mnt), ., value = T) %>% 
    as.character() %>% 
    rast() %>% 
    terra::crop(., zone) %>% 
    terra::mask(., zone) %>% 
    terra::crop(., c(-82, -66.9, -18.4, 12.5)) 
  srd <- srd/1000000
  
  if(nlyr(srd) == 29){
    print('Year leap')
    srd <- srd[[-29]]
  } else {
    print('Not leap')
    srd <- srd
  }
  
  # Temperature
  fls.tmx <- grep(glue('_{yr}'), grep('tmax', fls, value = T), value = T)
  fls.tmx <- grep(glue('_{yr}-{mnt}'), fls.tmx, value = T) 
  tmx <- rast(fls.tmx)
  
  fls.tmn <- grep(glue('_{yr}'), grep('tmin', fls, value = T), value = T)
  fls.tmn <- grep(glue('_{yr}-{mnt}'), fls.tmn, value = T) 
  tmn <- rast(fls.tmn)
  
  # Precipitation
  fls.ppt <- grep(glue('_{yr}'), grep('prec', fls, value = T), value = T)
  fls.ppt <- grep(glue('_{yr}-{mnt}'), fls.ppt, value = T)
  ppt <- rast(fls.ppt)
  
  # Raster reference
  ref <- ppt[[1]] * 0
  names(ref) <- 'ref'
  
  # To make the resample 
  tmx <- terra::resample(tmx, ref, method = 'bilinear')
  tmn <- terra::resample(tmn, ref, method = 'bilinear')
  srd <- terra::resample(srd, ref, method = 'bilinear')
  
  # To calculate average temperature
  tav <- (tmx + tmn) / 2
  
  # To calculate ETP
  ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
  
  # Compute water balance model
  date <- paste0(yr,'-',mnt)
  if(date %in% c('1986-01','2026-01')){
    cat('Start\n')
    AVAIL <<- ref
    AVAIL[!is.na(AVAIL)] <- 0
  } else {
    cat('Continous\n')
    AVAIL <<- terra::rast(glue('../data/tif/soils/avail.tif'))
  }
  
  # Fill gaps SCP - SST
  scp <- terra::focal(scp, w = 9, fun = mean, na.policy = 'only', na.rm=TRUE) 
  sst <- terra::focal(sst, w = 9, fun = mean, na.policy = 'only', na.rm=TRUE) 
  scp <- terra::crop(scp, zone) %>% terra::mask(., zone)
  sst <- terra::crop(sst, zone) %>% terra::mask(., zone)
  
  scp <- terra::resample(scp, ref)
  sst <- terra::resample(sst, ref)
  ppt <- terra::resample(ppt, ref)
  
  watbal <- map(.x = 1:nlyr(ETMAX), .f = function(i){
    wtr <- eabyep_calc(soilcp  = scp,
                       soilsat = sst,
                       avail   = AVAIL[[terra::nlyr(AVAIL)]],
                       rain    = ppt[[i]],
                       evap    = ETMAX[[i]])
    AVAIL <<- wtr$Availability
    return(wtr)
  })
  
  LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
  NDWL0  <- sum(LOGGING > 0) 
  names(NDWL0) <- glue('NDWL0_{yr}-{mn}_{md}')
  terra::writeRaster(x = AVAIL[[nlyr(AVAIL)]], glue('../data/tif/soils/avail.tif'), overwrite = TRUE)
  rm(ETMAX, LOGGING, ppt, tmx, tmn, tav, watbal)
  return(NDWL0)
  
}

# Load data ---------------------------------------------------------------

# Spatial data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Climate data 
root <- dir_ls('../data/tif/climate/future')
mdls <- dir_ls(root[1]) %>% dir_ls() %>% basename()

fles.tmax <- dir_ls('../data/tif/climate/baseline/tmax/daily', regexp = '.tif$')
fles.tmin <- dir_ls('../data/tif/climate/baseline/tmin/daily', regexp = '.tif$')
fles.prec <- dir_ls('../data/tif/climate/baseline/prec/daily', regexp = '.tif$')

# Solar radiation - Periods -----------------------------------------------
fles.srad <- dir_ls('../data/tif/srad/agERA5', type = 'directory') %>% as.character()
prds <- tibble(base = 1986:2015, future = 2026:2055)

# Soils -------------------------------------------------------------------
scp <- rast('../data/tif/soils/col_ecu_per_scp.tif') %>% terra::crop(., zone) %>% terra::mask(., zone)
sst <- rast('../data/tif/soils/col_ecu_per_ssat.tif')%>% terra::crop(., zone) %>% terra::mask(., zone)

# To calculate the NDWL0 --------------------------------------------------

years <- 2026:2055
month <- 1:12

# By one ------------------------------------------------------------------

# Model 1
ndwl.mdl1 <- map(.x = 1:length(years), .f = function(i){
  yea <- years[i]
  fnl <- map(.x = 1:length(month), .f = function(j){
    ndwl <- calc.ndwl0(md = mdls[1], yr = yea, mn = month[j])  
    return(ndwl)
  }) %>% 
    reduce(., c) %>% 
    app(., sum, na.rm = TRUE)
  names(fnl) <- glue('NDWL0_{yea}_{mdls[1]}')
  return(fnl)
}) %>% 
  reduce(., c)

# To write time series
dir_create('../data/tif/index/future/ndwl0')
terra::writeRaster(x = ndwl.mdl1, filename = '../data/tif/index/future/ndwl0/ndwl0_bsl_ACCESS-ESM1-5.tif', overwrite = TRUE)

# All into only one process -----------------------------------------------
ccl.mdl <- function(mdl){
  
  ndwl <- map(.x = 1:length(years), .f = function(i){
    yea <- years[i]
    fnl <- map(.x = 1:length(month), .f = function(j){
      ndwl <- calc.ndwl0(md = mdl, yr = yea, mn = month[j])  
      return(ndwl)
    }) %>% 
      reduce(., c) %>% 
      app(., sum, na.rm = TRUE)
    names(fnl) <- glue('NDWL0_{yea}_{mdl}')
    return(fnl)
  }) %>% 
    reduce(., c)
  
  terra::writeRaster(x = ndwl, filename = glue('../data/tif/index/future/ndwl0/ndwl0_bsl_{mdl}.tif'))
  cat('Finish!\n')
  
}
map(mdls[2:5], ccl.mdl)

