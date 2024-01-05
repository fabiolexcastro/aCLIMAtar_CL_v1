

# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# =========================================================================
# Function to use ---------------------------------------------------------
# =========================================================================
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.hsh.bsl <- function(yr){
  
  # yr <- 1990
  
  cat("To process: ", yr, '\n')
  
  # Files 
  fl.hm <- grep(paste0('_', yr), fles.hmdt, value = T)
  fl.tx <- grep(paste0('_', yr), fles.tmax, value = T)
  fl.tn <- grep(paste0('_', yr), fles.tmin, value = T)
  
  # As raster fiels and extract by mask 
  rs.hm <- rast(fl.hm) %>% crop(., zone) %>% mask(., zone)
  rs.tx <- rast(fl.tx) %>% crop(., zone) %>% mask(., zone)
  rs.tn <- rast(fl.tn) %>% crop(., zone) %>% mask(., zone)
  rs.hm <- terra::resample(rs.hm, rs.tx, method = 'bilinear')
  
  # rs.tx <- rs.tx - 273.15
  # rs.tn <- rs.tn - 273.15
  
  # Constant
  c1 = -8.78469475556
  c2 =  1.61139411
  c3 =  2.33854883889
  c4 = -0.14611605
  c5 = -0.012308094
  c6 = -0.0164248277778
  c7 =  2.211732 * 10^(-3)
  c8 =  7.2546 * 10^(-4)
  c9 = -3.582 * 10^(-6)
  
  # To calculate the average temperature 
  rs.ta <- (rs.tx + rs.tn) / 2
  
  hs <- purrr::map(.x = 1:12, .f = function(i){
    heat_idx <- function(ta, rh){
      hi <- ifelse(ta >= 25, c1 + (c2*ta) + (c3*rh) + (c4*ta*rh) + (c5*ta^2) + (c6*rh^2) + (c7*ta^2*rh) + (c8*ta*rh^2) + (c9*ta^2*rh^2), ta)
      return(hi)
    }
    hm <- rs.hm[[i]]
    ta <- rs.ta[[i]]
    HI <- terra::lapp(terra::sds(ta, hm), fun = heat_idx)
    HI_avg <- mean(HI, na.rm = T) %>% terra::mask(., zone)
    HI_max <- max(HI, na.rm = T) %>% terra::mask(., zone)
    return(list(HI_avg, HI_max))
  })
  
  hs.av <- do.call('c', map(hs, 1))
  hs.mx <- do.call('c', map(hs, 2))
  names(hs.av) <- glue('hs.avg_{yr}_{1:12}')
  names(hs.mx) <- glue('hs.max_{yr}_{1:12}')
  
  # plot(hs.av)
  cat('Done everything!\n')
  return(hs.av)
  
}

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
path <- '../data/tif/climate/baseline'
dirs.tasm <- dir_ls(path, regexp = 'tm') %>% as.character()

# Maximum temperature / Minimum temperature
fles.tmax <- grep('tmax', dirs.tasm, value = T) %>% dir_ls() %>% grep('monthly', ., value = T) %>% dir_ls() %>% as.character()
fles.tmin <- grep('tmin', dirs.tasm, value = T) %>% dir_ls() %>% grep('monthly', ., value = T) %>% dir_ls() %>% as.character()

# Humidity relative
dirs.hmdt <- dir_ls(path, regexp = 'rh') %>% as.character()
fles.hmdt <- dir_ls(dirs.hmdt) %>% grep('monthly', ., value = T) %>% dir_ls() %>% as.character()

# Spatial data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# =========================================================================
# To apply the function ---------------------------------------------------
# =========================================================================

# Time series
hsh.tsr <- map(.x = 1990:2015, .f = calc.hsh.bsl)
hsh.tsr <- reduce(hsh.tsr, c)

dir <- '../data/tif/index/baseline/hsh/hsh_bsl_tsr_chirts.tif'
dir_create(dirname(dir))
terra::writeRaster(x = hsh.tsr, filename = dir, overwrite = TRUE)




