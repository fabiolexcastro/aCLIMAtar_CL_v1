

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, future, furrr, lubridate, forecast, zoo, tidyverse, glue, stringr, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function to use ---------------------------------------------------------
make.tble <- function(var){
  
  # var <- 'tmax'
  cat('... Processing: ', var, '\n')
  
  fls <- dir_ls(path) %>% 
    grep(var, ., value = T) %>% 
    as.character() %>% 
    dir_ls(., regexp = 'monthly') %>% 
    dir_ls() %>% 
    as.character()
  
  rst <- map(fls, rast)
  names(rst)
  
  # Raster to table 
  tbl <- map(.x = 1:length(rst), .f = function(i){
    
    cat('To process: ', i, '\n')
    rs <- rst[[i]]
    fl <- fls[i]
    # nm <- basename(fl)
    # nm <- gsub('.tif', '', nm)
    # nm <- str_split(nm, pattern = '_') %>% map_chr(., 2)
    # nm <- paste0(var, '_', nm)
    # names(rs) <- nm
    
    tb <- terra::as.data.frame(rs, xy = T)
    tb <- as_tibble(tb)
    tb <- mutate(tb, gid = 1:nrow(tb))
    tb <- gather(tb, variable, value, -x, -y, -gid)
    tb <- mutate(tb, date = str_sub(variable, 6, nchar(variable)))
    tb <- mutate(tb, year = str_sub(date, 1, 4), month = str_sub(date, 6, nchar(date)))
    tb <- mutate(tb, year = as.numeric(year), month = as.numeric(month))
    tb <- mutate(tb, variable = str_sub(variable, 1, 4))
    tb <- dplyr::select(tb, gid, x, y, variable, year, month, value)
    return(tb)
    
  })
  
  tbl <- bind_rows(tbl)
  write.csv(tbl, glue('../data/tbl/historic/{var}_monthly.csv'), row.names = FALSE)
  cat('Done!\n')
  
}
make.frcs.1 <- function(tbl, iso, lbl){
  
  tbl <- tmin.tble
  iso <- 'COL'
  lbl <- 'tmin'
  
  cat('To start the process!\n')
  
  # To extract by mask
  dts <- tbl %>% distinct(year, month)
  lmt <- zone[zone$sov_a3 == iso,]
  lmt <- vect(lmt)
  
  trr <- map(.x = 1:nrow(dts), .f = function(d){
    
    cat('Raster: ', d, '\n')
    dt <- dts[d,]
    yr <- pull(dt, 1)
    mn <- pull(dt, 2)
    
    rs <- tbl %>% filter(year == yr & month == mn) %>% dplyr::select(x, y, value) %>% rast(., type = 'xyz') %>% terra::crop(., lmt) %>% terra::mask(., lmt)
    names(rs) <- glue('rs_{yr}-{mn}')
    return(rs)
    
  })
  
  trr <- reduce(trr, c)
  plot(trr[[1]])
  names(trr)
  terra::writeRaster(x = trr, filename = glue('../data/tif/climate/baseline/{lbl}/{lbl}_{iso}.tif'), overwrite = T)
  
  # Raster to table 
  dfm <- trr %>% 
    terra::as.data.frame(., xy = T) %>% 
    as_tibble() %>% 
    mutate(gid = 1:nrow(.)) %>% 
    gather(var, value, -x, -y, -gid) %>% 
    mutate(date = str_sub(var, 4, nchar(var))) %>% 
    dplyr::select(-var) %>% 
    separate(col = 'date', into = c('year', 'month'), sep = '-')
  
  gds <- dfm %>% 
    pull(gid) %>% 
    unique()
  
  options(future.globals.maxSize = 8000 * 1024^2)
  plan(cluster, workers = 18, gc = TRUE)
  rsl <- furrr::future_map(.x = gds, .f = function(i){
    
    cat('Pixel number ', i, '\n')
    tb <- filter(tbl, gid == i)
    
    fr <- map_dfr(.x = 1:12, .f = function(m){
      
      df <- tb %>% filter(month == m)
      vo <- unique(df$value)
      
      if(length(vo) == 1){
        cat('Unique value: ', 0, '\t')
        rs <- tibble(gid = i, month = m, date = glue('2023-{m}'), value = 0)
      } else {
        cat('Returned value not equal to 0\t')
      }
      
      vl <- pull(df, 7)
      tt <- ts(vl)
      ft <- auto.arima(tt)
      fc <- forecast(ft, h = 3)
      rs.mn <- fc$lower[,1]
      rs.mn <- as.numeric(rs.mn)
      rs.mx <- fc$upper[,1]
      rs.mx <- as.numeric(rs.mx)
      rs <- (rs.mn + rs.mx) / 2
      rs <- tibble(gid = i, month = m, date = c(glue('2023-{m}'), glue('2024-{m}'), glue('2025-{m}')), value = rs)
      # plot(c(vl, rs$value), type = 'l')
      return(rs)
      
    })
    
    # fr <- fr %>% group_by(gid, month) %>% dplyr::summarise(value = mean(value)) %>% ungroup()
    return(fr)
    
  })
  
  # To join all the individual tables into only one
  rsl2 <- rsl
  rsl2 <- bind_rows(rsl2)
  crds <- dfm %>% distinct(gid, x, y)
  
  # To check the number of rows
  pull(rsl2, gid) %>% unique() %>% length(); nrow(crds)
  
  # To join with the coordinates and saving as a rds file
  rsl2 <- inner_join(rsl2, crds, by = 'gid')
  
  if(lbl == 'prec'){
    print('Prec, checking < 0')
    rsl2 <- mutate(rsl2, value = ifelse(value < 0, 0, value))  
  }
  
  rsl2 %>% filter(date == '2023-1') %>% dplyr::select(x, y, value) %>% rast() %>% plot()
  
  saveRDS(object = rsl2, file = glue('../data/rds/forecast_{lbl}-{iso}.rds'))
  cat('Done!\n')
  
}
make.frcs.2 <- function(tbl, iso, lbl){
  
  tbl <- prec.tble  
  iso <- 'COL'
  lbl <- 'prec'
  
  terra::writeRaster(x = trr, filename = glue('../data/tif/climate/baseline/{lbl}/{lbl}_{iso}.tif'), overwrite = TRUE)
  
  rst <- dir_ls('../data/tif/climate/baseline') %>% 
    grep(lbl, ., value = T) %>% 
    as.character() %>% 
    dir_ls(., regexp = '.tif$') %>% 
    grep(iso, ., value = T) %>% 
    rast()
  
  # To check january
  rst.1 <- rst[[grep('-1$', names(rst), value = F)]]
  
  

}


# Load data ---------------------------------------------------------------

# Path raster data
path <- '../data/tif/climate/baseline'
vars <- c('prec', 'tmin', 'tmax')

# Administrative data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
isos <- c('COL', 'ECU', 'PER')
zone <- wrld[wrld$sov_a3 %in% isos,]

# Create the table --------------------------------------------------------

# Make table
prec.tble <- make.tble(var = 'prec')
tmax.tble <- make.tble(var = 'tmax')
tmin.tble <- make.tble(var = 'tmin')

# Read the result of the table
prec.tble <- read_csv('../data/tbl/historic/prec_monthly.csv')
tmax.tble <- read_csv('../data/tbl/historic/tmax_monthly.csv')
tmin.tble <- read_csv('../data/tbl/historic/tmin_monthly.csv')

# To apply the function ---------------------------------------------------

# 1. Make forescast processing - Way one

map(.x = 1:length(isos), .f = function(z){
  
  iso <- isos[z]
  
  cat('Processing: ', iso, '\n')
  make.frcs(tbl = prec.tble, iso = isos[z], lbl = 'prec')
  make.frcs(tbl = tmin.tble, iso = isos[z], lbl = 'tmin')
  make.frcs(tbl = tmax.tble, iso = isos[z], lbl = 'tmax')
  
})


# 2. Make forecast processing - Way two



# To check the results  ---------------------------------------------------





