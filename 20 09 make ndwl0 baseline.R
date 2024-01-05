
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

calc.ndwl0 <- function(yr, mn){
  
  # yr <- 1990
  # mn <- 1
  
  cat('To process: ', yr, ' ', mn, '\n')
  mnt <- mnth(mn)
  
  # Solar radiation
  srd <- rast(as.character(grep(glue('_{yr}{mnt}'), fles.srad, value = T)))
  srd <- terra::crop(srd, zone)
  srd <- terra::mask(srd, zone)
  srd <- srd/1000000
  
  # Temperature
  tmx <- rast(grep(glue('_{yr}'), fles.tmax, value = T))
  tmx <- tmx[[grep(glue('.{yr}.{mnt}.'), names(tmx), value = F)]]
  
  tmn <- rast(grep(glue('_{yr}'), fles.tmin, value = T))
  tmn <- tmn[[grep(glue('.{yr}.{mnt}.'), names(tmn), value = F)]]
  
  # tmx <- tmx - 273.15
  # tmn <- tmn - 273.15
  
  # Precipitation
  ppt <- rast(grep(glue('_{yr}'), fles.prec, value = T))
  ppt <- ppt[[grep(glue('{yr}.{mnt}'), names(ppt), value = F)]]
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
  if(date %in% c('1990-01','2026')){
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
fles.tmax <- dir_ls('../data/tif/climate/baseline/tmax/daily', regexp = '.tif$')
fles.tmin <- dir_ls('../data/tif/climate/baseline/tmin/daily', regexp = '.tif$')
fles.prec <- dir_ls('../data/tif/climate/baseline/prec/daily', regexp = '.tif$')
fles.srad <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/ERA5/solar_radiation_flux', regexp = '.nc$')

# Soils -------------------------------------------------------------------
scp <- rast('../data/tif/soils/col_ecu_per_scp.tif') %>% terra::crop(., zone) %>% terra::mask(., zone)
sst <- rast('../data/tif/soils/col_ecu_per_ssat.tif')%>% terra::crop(., zone) %>% terra::mask(., zone)

# To calculate the NDWL0 --------------------------------------------------

years <- 1990:2014
month <- 1:12

ndwl <- map(.x = 1:length(years), .f = function(i){
  yea <- years[i]
  fnl <- map(.x = 1:length(month), .f = function(j){
    ndwl <- calc.ndwl0(yr = yea, mn = month[j])  
    return(ndwl)
  }) %>% 
    reduce(., c) %>% 
    app(., sum, na.rm = TRUE)
  names(fnl) <- glue('NDWL0_{yea}')
  return(fnl)
}) %>% 
  reduce(., c)

# To write time series
terra::writeRaster(x = ndwl, filename = '../data/tif/indices/waf/ndwl0/ndwl0_bsl_tsr.tif')

ndwl <- terra::rast('../data/tif/indices/eaf/ndwl0/ndwl0_bsl_tsr.tif')



