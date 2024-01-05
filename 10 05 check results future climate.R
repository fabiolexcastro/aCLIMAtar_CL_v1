
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue)

rm(list = ls())
g <- gc(reset = T)
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
root <- dir_ls('../data/tif/climate/future', type = 'directory') %>% grep('/tm', ., value = T) %>% as.character()
sspe <- 'ssp370'
mdls <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
lbls <- rbind(tibble(year = 2021:2040, period = '2021-2040'), tibble(year = 2041:2060, period = '2041-2060'))

# Minimum temperature -----------------------------------------------------
tmin <- grep('tmin', root, value = T) %>% dir_ls() %>% dir_ls() %>% map(., dir_ls) %>% unlist() %>% as.character()
tmin.day <- grep('daily', tmin, value = T) %>% dir_ls(regexp = '.tif$') %>% as.character()

# Maximum temperature -----------------------------------------------------
tmax <- grep('tmax', root, value = T) %>% dir_ls() %>% dir_ls() %>% map(., dir_ls) %>% unlist() %>% as.character()
tmax.day <- grep('daily', tmax, value = T) %>% dir_ls(regexp = '.tif$') %>% as.character()

# Precipitation -----------------------------------------------------------
prec <- dir_ls('../data/tif/climate/future/prec') %>% dir_ls() %>% map(., dir_ls) %>% unlist() %>% as.character()
prec.day <- grep('daily', prec, value = T) %>% dir_ls(., regexp = '.tif$') %>% as.character()

# Function to use ---------------------------------------------------------
count.files <- function(fles){
  
  # To start the process
  print('Start the process')
  path <-  dirname(fles)
  
  # To make the table
  tble <- tibble(model = map_chr(str_split(fles, '/'), 8), file = basename(fles))
  tble <- separate(data = tble, col = 'file', into = c('variable', 'date'), sep = '_')
  tble <- mutate(tble, date = gsub('.tif$', '', date))
  tble <- separate(data = tble, col = 'date', into = c('year', 'month'), sep = '-')
  tble <- mutate(tble, year = as.numeric(year), month = as.numeric(month))
  tble <- inner_join(tble, lbls, by = 'year')
  
  # Get the summary for each period
  smmr.prdo <- tble %>% 
    group_by(model, variable, period) %>% 
    dplyr::summarise(count = n()) %>% 
    ungroup()
  
  smmr <- tble %>% 
    group_by(model, variable, year) %>% 
    dplyr::summarise(count = n()) %>% 
    ungroup()
  
  # tble %>% filter(model == 'ACCESS-ESM1-5') %>% group_by(model, period, year) %>% dplyr::summarise(count = n()) %>% ungroup() %>% View()

  print('Done!')
  return(smmr.prdo)
  
}

# To apply the function ---------------------------------------------------
smmr.tmin <- count.files(fles = tmin.day)
smmr.tmax <- count.files(fles = tmax.day)
smmr.prec <- count.files(fles = prec.day)
