

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, cowplot,  sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions ---------------------------------------------------------------
make.map.clm <- function(tbl, iso){
  
  # tbl <- tbls[[1]]
  # iso <- 'col'
  
  # Start the analysis 
  cat('To process:', iso, '\t')
  tbl <- mutate(tbl, gid = 1:nrow(tbl))
  tbl <- dplyr::select(tbl, 1:278)
  zne <- filter(zone, sov_a3 == toupper(iso))
  
  # Tidy the table ----------------------------------------------------------
  tbl <- tbl %>% 
    mutate(gid = 1:nrow(.)) %>% 
    gather(variable, value, -gid, -x, -y) %>% 
    separate(data = ., col = 'variable', sep = '_', into = c('variable', 'period', 'month')) %>% 
    filter(period %in% c('bsl', 'frc', 'ftr-avg', 'ftr')) %>% 
    mutate(period = factor(period, levels = c('bsl', 'frc', 'ftr-avg', 'ftr')), 
           month = as.numeric(month)) %>%
    inner_join(., tibble(month = 1:12, month_abb = month.abb), by = 'month') %>% 
    mutate(month_abb = factor(month_abb, levels = month.abb))
  
  unique(tbl$variable)
  
  # Maps for precipitation --------------------------------------------------
  tbl.ppt <- filter(tbl, variable == 'prec')
  
  gprec.bsl <- ggplot() +
    geom_tile(data = tbl.ppt %>% filter(period == 'bsl'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Precipitación - Línea base') +
    labs(fill = 'Prec (mm)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gprec.frc <- ggplot() +
    geom_tile(data = tbl.ppt %>% filter(period == 'frc'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    ggtitle(label = 'Precipitación - Actual') +
    labs(fill = 'Prec (mm)') +
    coord_sf() +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gprec.ftr <- ggplot() +
    geom_tile(data = tbl.ppt %>% filter(period == 'ftr-avg'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Precipitación - Futuro') +
    labs(fill = 'Prec (mm)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gprec <- ggarrange(gprec.bsl, gprec.frc, gprec.ftr, ncol = 3, nrow = 1, common.legend = TRUE)
  ggsave(plot = gprec, filename = glue('../png/maps/prec_{iso}.png'), units = 'in', width = 15, height = 7, dpi = 300)
  
  # Maps for average temperature --------------------------------------------
  tbl.tas <- filter(tbl, variable == 'tavg')
  
  gtasm.bsl <- ggplot() +
    geom_tile(data = tbl.tas %>% filter(period == 'bsl'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Temperatura promedio - Línea base') +
    labs(fill = 'Temperatura (°C)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gtasm.frc <- ggplot() +
    geom_tile(data = tbl.tas %>% filter(period == 'frc'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Temperatura promedio - Actual') +
    labs(fill = 'Temperatura (°C)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gtasm.ftr <- ggplot() +
    geom_tile(data = tbl.tas %>% filter(period == 'ftr'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Temperatura promedio - Futuro') +
    labs(fill = 'Temperatura (°C)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gtasm <- ggarrange(gtasm.bsl, gtasm.frc, gtasm.ftr, ncol = 3, nrow = 1, common.legend = TRUE)
  ggsave(plot = gtasm, filename = glue('../png/maps/tasm_{iso}.png'), units = 'in', width = 15, height = 7, dpi = 300)
  
  # Maps for ETP ------------------------------------------------------------
  tbl.etp <- filter(tbl, variable == 'etps')
  
  getp.bsl <- ggplot() +
    geom_tile(data = tbl.etp %>% filter(period == 'bsl'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = rev(brewer.pal(n = 9, name = 'BrBG'))) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    coord_sf() +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    ggtitle(label = 'Evapot. Potencial - Línea base') +
    labs(fill = 'ETP (mm)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  getp.frc <- ggplot() +
    geom_tile(data = tbl.etp %>% filter(period == 'frc'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = rev(brewer.pal(n = 9, name = 'BrBG'))) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Evapot. Potencial - Actual') +
    labs(fill = 'ETP (mm)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  getp.ftr <- ggplot() +
    geom_tile(data = tbl.etp %>% filter(period == 'ftr-avg'), aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = rev(brewer.pal(n = 9, name = 'BrBG'))) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Evapot. Potencial - Futuro') +
    labs(fill = 'ETP (mm)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          legend.key.width = unit(3, 'line'), 
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  getp <- ggarrange(getp.bsl, getp.frc, getp.ftr, ncol = 3, nrow = 1, common.legend = TRUE)
  ggsave(plot = getp, filename = glue('../png/maps/etps_{iso}.png'), units = 'in', width = 15, height = 7, dpi = 300)
  
  # Maps for balance --------------------------------------------------------
  tbl.bal <- filter(tbl, variable == 'baln')
  tbl.bal <- mutate(tbl.bal, value_cat = ifelse(value == 0, 'Normal', ifelse(value == -1, 'Dry', 'Wet')))
  tbl.bal <- mutate(tbl.bal, value_cat = factor(value_cat, levels = c('Dry', 'Normal', 'Wet')))
  
  gbal.bsl <- ggplot() +
    geom_tile(data = tbl.bal %>% filter(period == 'bsl'), aes(x = x, y = y, fill = value_cat)) + 
    scale_fill_manual(values = c('#99990C', '#F2F2F2', '#1D5C8F')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    coord_sf() +
    ggtitle(label = 'Balance - Línea base') +
    labs(fill = 'Balance (clase)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gbal.frc <- ggplot() +
    geom_tile(data = tbl.bal %>% filter(period == 'frc'), aes(x = x, y = y, fill = value_cat)) + 
    scale_fill_manual(values = c('#99990C', '#F2F2F2', '#1D5C8F')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    coord_sf() +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    ggtitle(label = 'Balance - Actual') +
    labs(fill = 'Balance (clase)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gbal.ftr <- ggplot() +
    geom_tile(data = tbl.bal %>% filter(period == 'ftr-avg'), aes(x = x, y = y, fill = value_cat)) + 
    scale_fill_manual(values = c('#99990C', '#F2F2F2', '#1D5C8F')) +
    facet_wrap(.~month_abb) + 
    # geom_sf(data = zne, fill = NA, col = 'grey30') +
    coord_sf() +
    # coord_sf(xlim = ext(zne)[1:2], ylim = ext(zne)[3:4]) +
    ggtitle(label = 'Balance - Futuro') +
    labs(fill = 'Balance (clase)') +
    theme_void() + 
    theme(legend.position = 'bottom', 
          plot.title = element_text(face = 'bold', hjust = 0.5),
          strip.text = element_text(face = 'bold'),
          plot.background = element_rect(colour = 'grey40', fill=NA, size = 1))
  
  gbal <- ggarrange(gbal.bsl, gbal.frc, gbal.ftr, ncol = 3, nrow = 1, common.legend = TRUE)
  ggsave(plot = gbal, filename = glue('../png/maps/baln_{iso}.png'), units = 'in', width = 15, height = 7, dpi = 300)
  
  cat('Finish!\n')
  
}
make.map.idx <- function(tbl, iso){
  
  # tbl <- tbls[[1]]
  # iso <- 'col'
  
  # Start the analysis 
  cat('>>>>> To process: ', iso, '<<<<<<\n')
  
  colnames(tbl)
  dfm <- tbl %>% 
    dplyr::select(
      x, y,
      ndd_bsl, ndd_ftr, tai_bsl, tai_ftr, ntx35_bsl, ntx35_ftr, hsh_bsl, hss_ftr, ndwl0_bsl, ndwl0_ftr
    ) %>% 
    mutate(
      gid = 1:nrow(.)
    ) %>% 
    gather(
      var, value, -c(gid, x, y)
    ) %>% 
    separate(
      col = 'var', into = c('index', 'period'), sep = '_'
    ) %>% 
    mutate(
      index = if_else(index == 'hss', 'hsh', index), 
      period = factor(period, levels = c('bsl', 'ftr'))
    ) 
  
  cls <- unique(dfm$value)
  
  # Number of dry days
  g.ndd <- ggplot() + 
    geom_tile(
      data = dfm %>% 
        filter(index == 'ndd') %>% 
        mutate(class = if_else(value == 1, 'Bajo', if_else(value == 2, 'Medio', 'Alto')),
               class = factor(class, levels = c('Bajo', 'Medio', 'Alto'))), 
      aes(
        x = x, 
        y = y, 
        fill = class
      )
    ) + 
    facet_wrap(
      .~period
    ) + 
    scale_fill_manual(
      values = brewer.pal(n = 3, name = 'YlOrBr')
    ) +
    coord_sf(
    ) +
    ggtitle(
      label = 'Numero de días secos'
    ) +
    theme_void(
    ) + 
    theme(
      legend.position = 'bottom', 
      strip.text = element_text(face = 'bold', hjust = 0.5, size = 14), 
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 17)
    ) + 
    labs(
      fill = 'Clase'
    )
  
  # Thornwaite Aridity Index
  g.tai <- ggplot() +
    geom_tile(
      data = dfm %>% 
        filter(index == 'tai') %>% 
        mutate(class = if_else(value == 1, 'Bajo', if_else(value == 2, 'Medio', 'Alto')),
               class = factor(class, levels = c('Bajo', 'Medio', 'Alto'))),
      aes(
        x = x,
        y = y, 
        fill = class
      )
    ) + 
    facet_wrap(
      .~period
    ) + 
    scale_fill_manual(
      values = brewer.pal(n = 3, name = 'YlOrBr')
    ) + 
    coord_sf(
    ) + 
    ggtitle(
      label = 'Índice de Aridez de Thornwaite'
    ) +
    theme_void(
    ) +
    theme(
      legend.position = 'bottom', 
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 16),
      strip.text = element_text(face = 'bold', hjust = 0.5, size = 14)
    ) + 
    labs(
      fill = 'Clases'
    )
  
  # NTX35
  g.ntx <- ggplot() + 
    geom_tile(
      data = dfm %>% 
        filter(index == 'ntx35') %>% 
        mutate(class = if_else(value == 1, 'Bajo', if_else(value == 2, 'Medio', 'Alto')),
               class = factor(class, levels = c('Bajo', 'Medio', 'Alto'))),
      aes(
        x = x, 
        y = y, 
        fill = class
      )
    ) +
    facet_wrap(
      .~period
    ) + 
    scale_fill_manual(
      values = brewer.pal(n = 3, name = 'YlOrRd')
    ) + 
    coord_sf(
    ) + 
    ggtitle(
      label = 'Número de días por encima de 35°C'
    ) +
    theme_void(
    ) + 
    theme(
      legend.position = 'bottom',
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 17),
      strip.text = element_text(face = 'bold', hjust = 0.5, size = 14)
    ) + 
    labs(
      fill = 'Clases'
    )
  
  # NDWL0
  g.ndw <- ggplot() + 
    geom_tile(
      data = dfm %>% 
        filter(index == 'ndwl0') %>% 
        mutate(class = if_else(value == 1, 'Bajo', if_else(value == 2, 'Medio', 'Alto')),
               class = factor(class, levels = c('Bajo', 'Medio', 'Alto'))),
      aes(
        x = x,
        y = y, 
        fill = class
      )
    ) +
    facet_wrap(
      .~period
    ) + 
    scale_fill_manual(
      values = brewer.pal(n = 3, name = 'YlOrRd')
    ) + 
    coord_sf(
    ) + 
    ggtitle(
      label = 'Número de días con el suelo inundado'
    ) +
    theme_void(
    ) + 
    theme(
      legend.position = 'bottom',
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 17),
      strip.text = element_text(face = 'bold', hjust = 0.5, size = 14)
    ) + 
    labs(
      fill = 'Clases'
    )
  
  # HSH
  g.hsh <- ggplot() +
    geom_tile(
      data = dfm %>% 
        filter(index == 'hsh') %>% 
        mutate(class = if_else(value == 1, 'Bajo', if_else(value == 2, 'Medio', if_else(value == 3, 'Alto', 'Muy alto'))),
               class = factor(class, levels = c('Bajo', 'Medio', 'Alto', 'Muy alto'))),
      aes(
        x = x,
        y = y, 
        fill = class
      )
    ) +
    facet_wrap(
      .~period
    ) + 
    scale_fill_manual(
      values = brewer.pal(n = 4, name = 'YlOrRd')
    ) + 
    coord_sf(
    ) + 
    ggtitle(
      label = 'Índice de calor humano'
    ) +
    theme_void(
    ) + 
    theme(
      legend.position = 'bottom',
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 17),
      strip.text = element_text(face = 'bold', hjust = 0.5, size = 14)
    ) + 
    labs(
      fill = 'Clases'
    )
  
  # Finish
  gall <- ggpubr::ggarrange(g.hsh, g.ndd, g.ndw, g.ntx, g.tai, nrow = 3, ncol = 2)
  ggsave(plot = gall, filename = glue('../png/maps/index_{iso}.png'), units = 'in', width = 11, height = 15, dpi = 300)
  cat('Done!\n')
  
  
}

# Load data ---------------------------------------------------------------
zone <- st_read('../data/gpk/zone.gpkg')
fles <- as.character(dir_ls('../data/tbl/results'))
tbls <- map(fles, read_csv)

# To make maps for the climate --------------------------------------------
map2(tbls, c('col', 'ecu', 'per'), make.map.clm)

# To make maps for the index ----------------------------------------------
map2(tbls, c('col', 'ecu', 'per'), make.map.idx)


