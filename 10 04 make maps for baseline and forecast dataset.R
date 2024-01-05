
## Fabio Castro - Llanos 
## Make maps for forecast results
## Alliance Bioversity - CIAT

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue, crayon, rnaturalearth)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
read.rds <- function(file){
  
  cat(green('Processing: ', basename(file), '\n'))
  
  # Read the rds file and make the summarise
  tble <- readRDS(file = file)
  varb <- basename(file) %>%
    str_split(., pattern = '_') %>% 
    map_chr(2) %>% 
    str_sub(., start = 1, end = 4)
  
  smmr <- tble %>% 
    group_by(gid, x, y, month) %>% 
    dplyr::summarise(value = mean(value, na.rm = T)) %>% 
    ungroup()
  
  
  # Conditional for check the negative values in precipitation variable
  if(varb == 'prec'){
    
    cat(yellow('Variable: precipitation\n'))
    smmr <- mutate(smmr, value = ifelse(value < 0, 0, value))
    
  } else { 
    
    cat(blue('Variable: temperature\n'))
    
  }
  
  # Convert table to raster
  smmr <- spread(data = smmr, key = month, value = value)
  rstr <- terra::rast(smmr[,2:ncol(smmr)])
  names(rstr) <- glue('{varb}_{1:12}')
  cat('Done!\n')
  return(rstr)
  
}


# Load data ---------------------------------------------------------------

# Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Tabular data 
fles <- dir_ls('../data/rds', regexp = '.rds$')
fles <- as.character(fles)

# Read the rds file and convert to raster  --------------------------------
prec <- fles %>% 
  grep('prec', ., value = T) %>% 
  map(.x = ., .f = read.rds)

plot(prec[[3]][[1]])


