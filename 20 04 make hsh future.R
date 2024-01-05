

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.hsh.ftr <- function(mdl){

  # List the files
  cat("To process: ", mdl, '\n')
  drs <- as.character(grep(mdl, dirs, value = T))
  yrs <- 2026:2055
  
  # To calculate the index
  hssa <- map(.x = 1:length(yrs), .f = function(y){
    
    cat('... Processing: ', yrs[y], '\n')
    
    # List tmax and tmin
    fl.tx <- grep('tmax', drs, value = TRUE) %>% 
      dir_ls(., regexp = 'monthly') %>% 
      dir_ls() %>% 
      grep(yrs[y], ., value = T) %>% 
      as.character()
    
    fl.tn <- grep('tmin', drs, value = TRUE) %>% 
      dir_ls(., regexp = 'monthly') %>% 
      dir_ls() %>% 
      grep(yrs[y], ., value = T) %>% 
      as.character()
    
    # Filter humidity relative
    yr.bs <- filter(prds, future == yrs[y])
    yr.ft <- pull(yr.bs, future)
    yr.bs <- pull(yr.bs, base)
    
    # Read as a raster file
    hr <- rast(grep(yr.bs, fl.hm, value = T))
    tx <- rast(fl.tx)
    tn <- rast(fl.tn)   
    
    # To calculate average temperature 
    ta <- (tx + tn) / 2
    
    # To resample humidity relative
    hr <- terra::resample(hr, ta)
    
    # Constants
    c1 = -8.78469475556
    c2 =  1.61139411
    c3 =  2.33854883889
    c4 = -0.14611605
    c5 = -0.012308094
    c6 = -0.0164248277778
    c7 =  2.211732 * 10^(-3)
    c8 =  7.2546 * 10^(-4)
    c9 = -3.582 * 10^(-6)
    
    # To calculate the heat index
    hs <- purrr::map(.x = 1:12, .f = function(i){
      
      cat('To process: ', month.abb[i], '\n')
      heat_idx <- function(ta, rh){
        hi <- ifelse(ta >= 25, c1 + (c2*ta) + (c3*rh) + (c4*ta*rh) + (c5*ta^2) + (c6*rh^2) + (c7*ta^2*rh) + (c8*ta*rh^2) + (c9*ta^2*rh^2), ta)
        return(hi)
      }
      
      hmd <- hr[[i]]
      tav <- ta[[i]]
      HI <- terra::lapp(terra::sds(tav, hmd), fun = heat_idx)
      HI_avg <- mean(HI, na.rm = T) %>% terra::mask(., zone)
      HI_max <- max(HI, na.rm = T) %>% terra::mask(., zone)
      return(list(HI_avg, HI_max))
    
    })
    
    # To calculate the index
    hs.av <- do.call('c', map(hs, 1))
    hs.mx <- do.call('c', map(hs, 2))
    names(hs.av) <- glue('hs.avg_{yrs[y]}_{1:12}')
    names(hs.mx) <- glue('hs.max_{yrs[y]}_{1:12}')
    
    # To return the files
    cat('Done everything!\n')
    return(hs.av)
    
  })
  hssa <- reduce(hssa, c)
  
  # To write the final raster (monthly)
  dout <- glue('../data/tif/index/future/hss/hss_{mdl}_monthly.tif')
  terra::writeRaster(x = hssa, filename = dout, overwrite = TRUE)
  rm(hssa); gc(reset = T)

}

# Load data ---------------------------------------------------------------
path <- '../data/tif/climate/future'
vars <- c('tmin', 'tmax')
dirs <- path %>% dir_ls() %>% dir_ls(., regexp = 'ssp370') %>% dir_ls() 

# Humidity files
fl.hm <- '../data/tif/climate/baseline/rhum' %>% dir_ls(., regexp = 'monthly') %>% dir_ls() %>% as.character()

# Make the period table
prds <- tibble(base = 1990:2019, future = 2026:2055)

# Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Models
mdls <- unique(basename(dirs))

# To calculate the index by each model ------------------------------------

# All in just one step
map(mdls, calc.hsh.ftr)

# Individually
calc.hsh.ftr(mdl = 'INM-CM5-0')
calc.hsh.ftr(mdl = 'MPI-ESM1-2-HR')
calc.hsh.ftr(mdl = 'MRI-ESM2-0')

# End ---------------------------------------------------------------------



