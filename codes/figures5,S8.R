rm(list = ls())

setwd("~/Documents/lab/body_size_evol")

library(phytools)
library(geiger)
library(RRphylo)
library(rgdal)
library(sf)
library(stringi)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(fasterize)

data(wrld_simpl, package = "maptools")

tr_am <- read.nexus("data/amphibia/MCC_amphibia.nex")
tr_sq <- read.nexus("data/reptilia/MCC_squamata.nex")
tr_av <- read.nexus("data/aves/MCC_aves.nex")
tr_ma <- read.nexus("data/mammalia/MCC_mammalia.nex")

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

mass_am <- dat_am$Body_mass_g_1[complete.cases(dat_am$Body_mass_g_1)]
names(mass_am) <- dat_am$Scientific_name[complete.cases(dat_am$Body_mass_g_1)]

svl_am <- dat_am$SVL_mm_1[complete.cases(dat_am$SVL_mm_1)]
names(svl_am) <- dat_am$Scientific_name[complete.cases(dat_am$SVL_mm_1)]

mass_sq <- dat_sq$Body_mass_g_mean[complete.cases(dat_sq$Body_mass_g_mean)]
names(mass_sq) <- dat_sq$Scientific_name[complete.cases(dat_sq$Body_mass_g_mean)]

svl_sq <- dat_sq$SVL_mm_mean[complete.cases(dat_sq$SVL_mm_mean)]
names(svl_sq) <- dat_sq$Scientific_name[complete.cases(dat_sq$SVL_mm_mean)]

mass_av <- dat_av$Body_mass_g_mean[complete.cases(dat_av$Body_mass_g_mean)]
names(mass_av) <- dat_av$Scientific_name[complete.cases(dat_av$Body_mass_g_mean)]

mass_ma <- dat_ma$Body_mass_g_mean[complete.cases(dat_ma$Body_mass_g_mean)]
names(mass_ma) <- dat_ma$Scientific_name[complete.cases(dat_ma$Body_mass_g_mean)]

pruned_tr_am <- treedata(tr_am, mass_am, warnings = F)$phy
pruned_tr_sq <- treedata(tr_sq, mass_sq, warnings = F)$phy
pruned_tr_av <- treedata(tr_av, mass_av, warnings = F)$phy
pruned_tr_ma <- treedata(tr_ma, mass_ma, warnings = F)$phy

pruned_tr2_am <- treedata(tr_am, svl_am, warnings = F)$phy
pruned_tr2_sq <- treedata(tr_sq, svl_sq, warnings = F)$phy

mass_am <- mass_am[names(mass_am) %in% pruned_tr_am$tip.label]
mass_sq <- mass_sq[names(mass_sq) %in% pruned_tr_sq$tip.label]
mass_av <- mass_av[names(mass_av) %in% pruned_tr_av$tip.label]
mass_ma <- mass_ma[names(mass_ma) %in% pruned_tr_ma$tip.label]

svl_am <- svl_am[names(svl_am) %in% pruned_tr2_am$tip.label]
svl_sq <- svl_sq[names(svl_sq) %in% pruned_tr2_sq$tip.label]

mass_am <- mass_am[match(pruned_tr_am$tip.label, names(mass_am))]
mass_sq <- mass_sq[match(pruned_tr_sq$tip.label, names(mass_sq))]
mass_av <- mass_av[match(pruned_tr_av$tip.label, names(mass_av))]
mass_ma <- mass_ma[match(pruned_tr_ma$tip.label, names(mass_ma))]

svl_am <- svl_am[match(pruned_tr2_am$tip.label, names(svl_am))]
svl_sq <- svl_sq[match(pruned_tr2_sq$tip.label, names(svl_sq))]

log_mass_am <- log(mass_am)
log_mass_sq <- log(mass_sq)
log_mass_av <- log(mass_av)
log_mass_ma <- log(mass_ma)

log_svl_am <- log(svl_am)
log_svl_sq <- log(svl_sq)

RR_am <- RRphylo(tree = pruned_tr_am, y = log_mass_am)
#save(RR_am, file = "data/amphibia/RR_am.RData")

RR_sq <- RRphylo(tree = pruned_tr_sq, y = log_mass_sq)
#save(RR_sq, file = "data/reptilia/RR_sq.RData")

RR_av <- RRphylo(tree = pruned_tr_av, y = log_mass_av)
#save(RR_av, file = "data/aves/RR_av.RData")

RR_ma <- RRphylo(tree = pruned_tr_ma, y = log_mass_ma)
#save(RR_ma, file = "data/mammalia/RR_ma.RData")

RR2_am <- RRphylo(tree = pruned_tr2_am, y = log_svl_am)
#save(RR2_am, file = "data/amphibia/RR2_am.RData")

RR2_sq <- RRphylo(tree = pruned_tr2_sq, y = log_svl_sq)
#save(RR2_sq, file = "data/reptilia/RR2_sq.RData")

#load(file = "data/amphibia/RR2_am.RData")
#load(file = "data/reptilia/RR2_sq.RData")

map_am <- readOGR("~/Documents/lab/data/spatial/AMPHIBIANS/AMPHIBIANS.shp",
                  dropNULLGeometries = TRUE)
map_sq1 <- readOGR("~/Documents/lab/data/spatial/SQUAMATA/SCALED_REPTILES_PART1.shp", 
                   dropNULLGeometries = TRUE)
map_sq2 <- readOGR("~/Documents/lab/data/spatial/SQUAMATA/SCALED_REPTILES_PART2.shp",
                   dropNULLGeometries = TRUE)
map_av <- st_read(dsn = "~/Documents/lab/data/spatial/BOTW/BOTW.gdb", 
                  layer = "All_Species")
map_ma <- readOGR("~/Documents/lab/data/spatial/MAMMALS/MAMMALS.shp",
                  dropNULLGeometries = TRUE)

map_av <- map_av %>% filter(st_geometry_type(Shape) != "MULTISURFACE")

map_am$sci_name <- stri_replace_all_fixed(map_am$sci_name, " ", "_")
map_sq1$sci_name <- stri_replace_all_fixed(map_sq1$sci_name, " ", "_")
map_sq2$sci_name <- stri_replace_all_fixed(map_sq2$sci_name, " ", "_")
map_av$sci_name <- stri_replace_all_fixed(map_av$sci_name, " ", "_")
map_ma$sci_name <- stri_replace_all_fixed(map_ma$sci_name, " ", "_")

make_raster <- function(maps, trait, taxon, path) {

  if(taxon == "Aves") {
    
    r <- raster(ncols = 2160, nrows = 900, ymn = -60)
    stack <- r
    
    for (i in 1:length(trait)) {
      print(paste0(i, "/", length(trait)))
      
      s <- as.character(names(trait)[i]) 
      map_i <- subset(maps, maps$sci_name == s)
      
      for (j in 1:length(map_i$Shape)) { 
        try(map_i$Shape[[j]] <- st_cast(map_i$Shape[[j]], 'MULTIPOLYGON'))
      }
      try(map_i$Shape <- st_cast(map_i$Shape, 'MULTIPOLYGON'))
      
      raster_i <- fasterize(st_as_sf(map_i$Shape), r)
      
      rastercells <- which(getValues(!is.na(raster_i))) 
      raster_i[rastercells] <- trait[i]
      
      writeRaster(raster_i, paste0(path, "/raster_", s, ".tif"), format = "GTiff")
    }
    
  } else {
    
    r <- raster(ncols = 2160, nrows = 900, ymn = -60)
    stack <- r
  
    for (i in 1:length(trait)) {
      print(paste0(i, "/", length(trait)))
    
      s <- as.character(names(trait)[i]) 
      map_i <- subset(maps, maps$sci_name == s)
    
      raster_i <- rasterize(map_i, r)
    
      rastercells <- which(getValues(!is.na(raster_i))) 
      raster_i[rastercells] <- trait[i]
    
      writeRaster(raster_i, paste0(path, "/raster_", s, ".tif"), format = "GTiff")
    }
  }
}

crop_rg <- function(path, extent, write, log = FALSE) {
 
  fList <- list.files(path, ".tif", full.names = T)

  for (i in 1:length(fList)) {
    print(paste0(i, "/", length(fList)))

    raster_i <- raster(fList[i])

    xx <- strsplit(names(raster_i), "_")[[1]]

    s <- paste0(xx[2], "_", xx[3])

    crop_i <- crop(raster_i, extent)
    
    if(log) {
      values(crop_i) <- log(values(crop_i))
    }
    if(!all(is.na(values(crop_i)))) {
      writeRaster(crop_i, paste0(write, "/raster_", s, ".tif"), 
                  format = "GTiff")
    }
    rm(raster_i)
    rm(crop_i)
    removeTmpFiles(h = 0)
  }
}

# Figure 5

rates_am <- RR_am$rates[, 1][names(RR_am$rates[, 1]) %in% tr_am$tip.label]
rates_sq <- RR_sq$rates[, 1][names(RR_sq$rates[, 1]) %in% tr_sq$tip.label]
rates_av <- RR_av$rates[, 1][names(RR_av$rates[, 1]) %in% tr_av$tip.label]
rates_ma <- RR_ma$rates[, 1][names(RR_ma$rates[, 1]) %in% tr_ma$tip.label]

map2_am <- raster::subset(map_am, map_am$sci_name %in% names(rates_am))
map2_sq1 <- raster::subset(map_sq1, map_sq1$sci_name %in% names(rates_sq))
map2_sq2 <- raster::subset(map_sq2, map_sq2$sci_name %in% names(rates_sq))
map2_av <- raster::subset(map_av, map_av$sci_name %in% names(rates_av))
map2_ma <- raster::subset(map_ma, map_ma$sci_name %in% names(rates_ma))

rates_am <- rates_am[names(rates_am) %in% map2_am$sci_name]
rates_sq1 <- rates_sq[names(rates_sq) %in% map2_sq1$sci_name]
rates_sq2 <- rates_sq[names(rates_sq) %in% map2_sq2$sci_name]
rates_av <- rates_av[names(rates_av) %in% map2_av$sci_name]
rates_ma <- rates_ma[names(rates_ma) %in% map2_ma$sci_name]

mass_am <- mass_am[names(mass_am) %in% map2_am$sci_name]
mass_sq1 <- mass_sq[names(mass_sq) %in% map2_sq1$sci_name]
mass_sq2 <- mass_sq[names(mass_sq) %in% map2_sq2$sci_name]
mass_av <- mass_av[names(mass_av) %in% map2_av$sci_name]
mass_ma <- mass_ma[names(mass_ma) %in% map2_ma$sci_name]

make_raster(map2_am, mass_am, "Amphibia", "/Volumes/Personal/lab/raster/amphibia_mass")
make_raster(map2_am, rates_am,"Amphibia", "/Volumes/Personal/lab/raster/amphibia_rates")

make_raster(map2_sq1, mass_sq1, "Squamata", "/Volumes/Personal/lab/raster/squamata_mass")
make_raster(map2_sq2, mass_sq2, "Squamata", "/Volumes/Personal/lab/raster/squamata_mass")
make_raster(map2_sq1, rates_sq1, "Squamata", "/Volumes/Personal/lab/raster/squamata_rates")
make_raster(map2_sq2, rates_sq2, "Squamata", "/Volumes/Personal/lab/raster/squamata_rates")

make_raster(map2_av, mass_av, "Aves", "/Volumes/Personal/lab/raster/aves_mass")
make_raster(map2_av, rates_av, "Aves", "/Volumes/Personal/lab/raster/aves_mass")

make_raster(map2_ma, mass_ma, "Mammalia", "/Volumes/Personal/lab/raster/mammalia_mass")
make_raster(map2_ma, rates_ma, "Mammalia", "/Volumes/Personal/lab/raster/mammalia_rates")

r_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
               xmn = -180, xmx = 180)
r1_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -180, xmx = -120)
r2_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -119.9999, xmx = -60)
r3_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -59.9999, xmx = 0.9999)
r4_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 1, xmx = 60.9999)
r5_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 61, xmx = 120.9999)
r6_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 121, xmx = 180)

writeRaster(r_am, "/Volumes/Personal/lab/raster/amphibia_mass/r.tif",
            format = "GTiff")

writeRaster(r1_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext1", log = F)

writeRaster(r2_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(-119.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext2", log = F)

writeRaster(r3_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(-59.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext3", log = F)

writeRaster(r4_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext4", log = F)

writeRaster(r5_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext5", log = F)

writeRaster(r6_am, "/Volumes/Personal/lab/raster/amphibia_mass/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_mass", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_mass/ext6", log = F)

writeRaster(r_am, "/Volumes/Personal/lab/raster/amphibia_rates/r.tif",
            format = "GTiff")

writeRaster(r1_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext1", log = F)

writeRaster(r2_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(-119.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext2", log = F)

writeRaster(r3_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(-59.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext3", log = F)

writeRaster(r4_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext4", log = F)

writeRaster(r5_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext5", log = F)

writeRaster(r6_am, "/Volumes/Personal/lab/raster/amphibia_rates/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates/ext6", log = F)

r_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
               xmn = -180, xmx = 180)
r1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -180, xmx = -120)
r2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -119.9999, xmx = -40)
r3_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -39.9999, xmx = 0.9999)
r4_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 1, xmx = 60.9999)
r5_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 61, xmx = 120.9999)
r6_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 121, xmx = 180)

writeRaster(r_sq, "/Volumes/Personal/lab/raster/squamata_mass/r.tif",
            format = "GTiff")

writeRaster(r1_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext1", log = F)

writeRaster(r2_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(-119.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext2", log = F)

writeRaster(r3_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext3", log = F)

writeRaster(r4_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext4", log = F)

writeRaster(r5_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext5", log = F)

writeRaster(r6_sq, "/Volumes/Personal/lab/raster/squamata_mass/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_mass", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_mass/ext6", log = F)

writeRaster(r_sq, "/Volumes/Personal/lab/raster/squamata_rates/r.tif",
            format = "GTiff")

writeRaster(r1_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext1", log = F)

writeRaster(r2_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(-119.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext2", log = F)

writeRaster(r3_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext3", log = F)

writeRaster(r4_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext4", log = F)

writeRaster(r5_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext5", log = F)

writeRaster(r6_sq, "/Volumes/Personal/lab/raster/squamata_rates/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates/ext6", log = F)


r_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
               xmn = -180, xmx = 180)
r1_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -180, xmx = -120)
r2_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -119.9999, xmx = -80)
r3_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -79.9999, xmx = -60)
r4_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -59.9999, xmx = -20)
r5_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -19.9999, xmx = 20.9999)
r6_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 21, xmx = 60.9999)
r7_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 61, xmx = 90.9999)
r8_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 91, xmx = 120.9999)
r9_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 121, xmx = 150.9999)
r10_av <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                 xmn = 151, xmx = 180)

writeRaster(r_av, "/Volumes/Personal/lab/raster/aves_mass/r.tif",
            format = "GTiff")

writeRaster(r1_av, "/Volumes/Personal/lab/raster/aves_mass/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext1", log = F)

writeRaster(r2_av, "/Volumes/Personal/lab/raster/aves_mass/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(-119.9999, -80, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext2", log = F)

writeRaster(r3_av, "/Volumes/Personal/lab/raster/aves_mass/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(-79.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext3", log = F)

writeRaster(r4_av, "/Volumes/Personal/lab/raster/aves_mass/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(-59.9999, -20, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext4", log = F)

writeRaster(r5_av, "/Volumes/Personal/lab/raster/aves_mass/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(-19.9999, 20.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext5", log = F)

writeRaster(r6_av, "/Volumes/Personal/lab/raster/aves_mass/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(21, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext6", log = F)

writeRaster(r7_av, "/Volumes/Personal/lab/raster/aves_mass/ext7/r7.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(61, 90.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext7", log = F)

writeRaster(r8_av, "/Volumes/Personal/lab/raster/aves_mass/ext8/r8.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(91, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext8", log = F)

writeRaster(r9_av, "/Volumes/Personal/lab/raster/aves_mass/ext9/r9.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(121, 150.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext9", log = F)

writeRaster(r10_av, "/Volumes/Personal/lab/raster/aves_mass/ext10/r10.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_mass", 
        extent(151, 180, -90, 90),
        "/Volumes/Personal/lab/raster/aves_mass/ext10", log = F)


writeRaster(r_av, "/Volumes/Personal/lab/raster/aves_rates/r.tif",
            format = "GTiff")

writeRaster(r1_av, "/Volumes/Personal/lab/raster/aves_rates/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext1", log = F)

writeRaster(r2_av, "/Volumes/Personal/lab/raster/aves_rates/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-119.9999, -80, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext2", log = F)

writeRaster(r3_av, "/Volumes/Personal/lab/raster/aves_rates/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-79.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext3", log = F)

writeRaster(r4_av, "/Volumes/Personal/lab/raster/aves_rates/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-59.9999, -20, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext4", log = F)

writeRaster(r5_av, "/Volumes/Personal/lab/raster/aves_rates/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-19.9999, 20.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext5", log = F)

writeRaster(r6_av, "/Volumes/Personal/lab/raster/aves_rates/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(21, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext6", log = F)

writeRaster(r7_av, "/Volumes/Personal/lab/raster/aves_rates/ext7/r7.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(61, 90.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext7", log = F)

writeRaster(r8_av, "/Volumes/Personal/lab/raster/aves_rates/ext8/r8.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(91, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext8", log = F)

writeRaster(r9_av, "/Volumes/Personal/lab/raster/aves_rates/ext9/r9.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(121, 150.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext9", log = F)

writeRaster(r10_av, "/Volumes/Personal/lab/raster/aves_rates/ext10/r10.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(151, 180, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext10", log = F)

writeRaster(r_av, "/Volumes/Personal/lab/raster/aves_rates/r.tif",
            format = "GTiff")

writeRaster(r1_av, "/Volumes/Personal/lab/raster/aves_rates/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext1", log = F)

writeRaster(r2_av, "/Volumes/Personal/lab/raster/aves_rates/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-119.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext2", log = F)

writeRaster(r3_av, "/Volumes/Personal/lab/raster/aves_rates/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext3", log = F)

writeRaster(r4_av, "/Volumes/Personal/lab/raster/aves_rates/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext4", log = F)

writeRaster(r5_av, "/Volumes/Personal/lab/raster/aves_rates/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext5", log = F)

writeRaster(r6_av, "/Volumes/Personal/lab/raster/aves_rates/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/aves_rates", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/aves_rates/ext6", log = F)



r_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = -180, xmx = 180)
r1_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = -180, xmx = -120)
r2_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = -119.9999, xmx = -40)
r3_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = -39.9999, xmx = 0.9999)
r4_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = 1, xmx = 60.9999)
r5_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = 61, xmx = 120.9999)
r6_ma <- raster(ncols = 2160, nrows = 900, ymn = -90, 
             xmn = 121, xmx = 180)

writeRaster(r_ma, "/Volumes/Personal/lab/raster/mammalia_mass/r.tif",
            format = "GTiff")

writeRaster(r1_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext1", log = F)

writeRaster(r2_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(-119.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext2", log = F)

writeRaster(r3_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext3", log = F)

writeRaster(r4_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext4", log = F)

writeRaster(r5_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext5", log = F)

writeRaster(r6_ma, "/Volumes/Personal/lab/raster/mammalia_mass/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_mass", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_mass/ext6", log = F)

writeRaster(r_ma, "/Volumes/Personal/lab/raster/mammalia_rates/r.tif",
            format = "GTiff")

writeRaster(r1_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext1", log = F)

writeRaster(r2_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(-119.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext2", log = F)

writeRaster(r3_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext3", log = F)

writeRaster(r4_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext4", log = F)

writeRaster(r5_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext5", log = F)

writeRaster(r6_ma, "/Volumes/Personal/lab/raster/mammalia_rates/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/mammalia_rates", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/mammalia_rates/ext6", log = F)

get_values <- function(full) {

  data(wrld_simpl, package = "maptools")

  rem1 <- extract(full, wrld_simpl, cellnumbers = T, weights = T, 
                  small = T)
  rem2 <- do.call(rbind.data.frame, rem1)[, 1]
  values(full)[-rem2] <- NA

  brks <- quantile(values(full)[order(values(full))], probs = seq(0, 1, 0.02),
                   na.rm = T)
  brks <- unique(brks)
  full_p <- rasterToPoints(full, spatial = TRUE)
  full_df  <- data.frame(full_p)
  colnames(full_df)[1] <- "index_1"
  full_df <- full_df %>% mutate(index_2 = cut(index_1, breaks = brks))

  return(full_df)
}

raster_mass_am <- get_values(raster("data/amphibia/mass.tif"))
raster_rates_am <- get_values(raster("data/amphibia/rates.tif"))

raster_mass_sq <- get_values(raster("data/reptilia/mass.tif"))
raster_rates_sq <- get_values(raster("data/reptilia/rates.tif"))

raster_mass_av <- get_values(raster("data/aves/mass.tif"))
raster_rates_av <- get_values(raster("data/aves/rates.tif"))

raster_mass_ma <- get_values(raster("data/mammalia/mass.tif"))
raster_rates_ma <- get_values(raster("data/mammalia/rates.tif"))

r <- raster(ncols = 2160, nrows = 900, ymn = -60)

colors_mass <- mako(50)[50:1]
colors_rates <- rocket(50)[50:1]

ras_wrld <- rasterToPoints(rasterize(wrld_simpl, r), spatial = TRUE)
ras_wrld_df  <- data.frame(ras_wrld)

plot1_am <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_mass_am, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mass) +
  theme_void() 


plot2_am <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates_am, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

plot1_sq <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_mass_sq, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mass) +
  theme_void() 


plot2_sq <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates_sq, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

plot1_av <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_mass_av, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mass) +
  theme_void() 


plot2_av <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates_av, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

plot1_ma <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_mass_ma, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mass) +
  theme_void() 


plot2_ma <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates_ma, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

inches <- 4.5
res <- 600

tiff("figures/Figure5_am1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_am)  
dev.off() 

tiff("figures/Figure5_am2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_am)  
dev.off() 

tiff("figures/Figure5_sq1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_sq)  
dev.off() 

tiff("figures/Figure5_sq2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_sq)  
dev.off() 

tiff("figures/Figure5_av1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_av)  
dev.off() 

tiff("figures/Figure5_av2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_av)  
dev.off() 

tiff("figures/Figure5_ma1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_ma)  
dev.off() 

tiff("figures/Figure5_ma2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_ma)  
dev.off() 

tiff("figures/Figure5_ScaleBar1.tiff", width = 1*res, h = 0.1*res, 
     units = "px")
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))
brks <- seq(0, 1, 0.02)
breaks <- seq(0, 100, length.out = length(brks))

ix <- 1:2
iy <- breaks
nBreaks <- length(breaks)
midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      col = colors_mass, breaks = breaks)
axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                  cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                  labels = rep("", 5), tck = 0.2)
do.call("axis", axis.args)
box()
dev.off()

tiff("figures/Figure5_ScaleBar2.tiff", width = 1*res, h = 0.1*res, 
     units = "px")
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))
brks <- seq(0, 1, 0.02)
breaks <- seq(0, 100, length.out = length(brks))

ix <- 1:2
iy <- breaks
nBreaks <- length(breaks)
midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      col = colors_rates, breaks = breaks)
axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                  cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                  labels = rep("", 5), tck = 0.2)
do.call("axis", axis.args)
box()
dev.off()

quantile(round(unique(quantile(values(raster("data/amphibia/mass.tif"))[order(values(raster("data/amphibia/mass.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/amphibia/rates.tif"))[order(values(raster("data/amphibia/rates.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/reptilia/mass.tif"))[order(values(raster("data/reptilia/mass.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/reptilia/rates.tif"))[order(values(raster("data/reptilia/rates.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/aves/mass.tif"))[order(values(raster("data/aves/mass.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/aves/rates.tif"))[order(values(raster("data/aves/rates.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/mammalia/mass.tif"))[order(values(raster("data/mammalia/mass.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/mammalia/rates.tif"))[order(values(raster("data/mammalia/rates.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

# Figure S8

rates2_am <- RR2_am$rates[, 1][names(RR2_am$rates[, 1]) %in% tr_am$tip.label]
rates2_sq <- RR2_sq$rates[, 1][names(RR2_sq$rates[, 1]) %in% tr_sq$tip.label]

map2_2_am <- raster::subset(map_am, map_am$sci_name %in% names(rates2_am))
map2_2_sq1 <- raster::subset(map_sq1, map_sq1$sci_name %in% names(rates2_sq))
map2_2_sq2 <- raster::subset(map_sq2, map_sq2$sci_name %in% names(rates2_sq))

rates2_am <- rates2_am[names(rates2_am) %in% map2_2_am$sci_name]
rates2_sq1 <- rates2_sq[names(rates2_sq) %in% map2_2_sq1$sci_name]
rates2_sq2 <- rates2_sq[names(rates2_sq) %in% map2_2_sq2$sci_name]

svl_am <- svl_am[names(svl_am) %in% map2_2_am$sci_name]
svl_sq1 <- svl_sq[names(svl_sq) %in% map2_2_sq1$sci_name]
svl_sq2 <- svl_sq[names(svl_sq) %in% map2_2_sq2$sci_name]

make_raster(map2_2_am, svl_am, "Amphibia", "/Volumes/Personal/lab/raster/amphibia_svl")
make_raster(map2_2_am, rates2_am,"Amphibia", "/Volumes/Personal/lab/raster/amphibia_rates_svl")

make_raster(map2_2_sq1, svl_sq1, "Squamata", "/Volumes/Personal/lab/raster/squamata_svl")
make_raster(map2_2_sq2, svl_sq2, "Squamata", "/Volumes/Personal/lab/raster/squamata_svl")

make_raster(map2_2_sq1, rates2_sq1[4196:length(rates2_sq1)], "Squamata", "/Volumes/Personal/lab/raster/squamata_rates_svl")
make_raster(map2_2_sq2, rates2_sq2, "Squamata", "/Volumes/Personal/lab/raster/squamata_rates_svl")

r_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
               xmn = -180, xmx = 180)
r1_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -180, xmx = -120)
r2_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -119.9999, xmx = -60)
r2.1_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -119.9999, xmx = -90)
r2.2_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -89.9999, xmx = -75)
r2.3_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -74.9999, xmx = -60)
r3_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -59.9999, xmx = 0.9999)
r4_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 1, xmx = 60.9999)
r5_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 61, xmx = 120.9999)
r6_am <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 121, xmx = 180)

writeRaster(r_am, "/Volumes/Personal/lab/raster/amphibia_svl/r.tif",
            format = "GTiff")

writeRaster(r1_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext1", log = F)

writeRaster(r2_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext2/r2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(-119.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext2", log = F)

writeRaster(r3_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(-59.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext3", log = F)

writeRaster(r4_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext4", log = F)

writeRaster(r5_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext5", log = F)

writeRaster(r6_am, "/Volumes/Personal/lab/raster/amphibia_svl/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_svl", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_svl/ext6", log = F)

writeRaster(r_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/r.tif",
            format = "GTiff")

writeRaster(r1_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext1", log = F)

writeRaster(r2.1_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.1/r2.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(-119.9999, -90, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.1", log = F)

writeRaster(r2.2_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.2/r2.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(-89.9999, -75, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.2", log = F)

writeRaster(r2.3_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.3/r2.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(-74.9999, -60, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext2.3", log = F)

writeRaster(r3_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(-59.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext3", log = F)

writeRaster(r4_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext4/r4.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(1, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext4", log = F)

writeRaster(r5_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext5/r5.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(61, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext5", log = F)

writeRaster(r6_am, "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext6/r6.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/amphibia_rates_svl", 
        extent(121, 180, -90, 90),
        "/Volumes/Personal/lab/raster/amphibia_rates_svl/ext6", log = F)

r_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
               xmn = -180, xmx = 180)
r1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -180, xmx = -120)
r2.1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -119.9999, xmx = -100)
r2.2.1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -99.9999, xmx = -85)
r2.2.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -84.9999, xmx = -75)
r2.3_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                  xmn = -74.9999, xmx = -40)
r3_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = -39.9999, xmx = 0.9999)
r4.1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 1, xmx = 30.9999)
r4.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 31, xmx = 45.9999)
r4.3_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 46, xmx = 60.9999)
r5.1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 61, xmx = 90.9999)
r5.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 91, xmx = 105.9999)
r5.3_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 106, xmx = 120.9999)
r6.1.1.1_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 121, xmx = 128.9999)
r6.1.1.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 129, xmx = 135.9999)
r6.1.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 136, xmx = 150.9999)
r6.2_sq <- raster(ncols = 2160, nrows = 900, ymn = -90, 
                xmn = 151, xmx = 180)

writeRaster(r_sq, "/Volumes/Personal/lab/raster/squamata_svl/r.tif",
            format = "GTiff")

writeRaster(r1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext1", log = F)

writeRaster(r2.1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext2.1/r2.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-119.9999, -100, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext2.1", log = F)

writeRaster(r2.2.1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext2.2.1/r2.2.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-99.9999, -85, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext2.2.1", log = F)

writeRaster(r2.2.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext2.2.2/r2.2.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-84.9999, -75, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext2.2.2", log = F)

writeRaster(r2.3_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext2.3/r2.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-74.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext2.3", log = F)

writeRaster(r3_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext3", log = F)

writeRaster(r4.1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext4.1/r4.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(1, 30.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext4.1", log = F)

writeRaster(r4.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext4.2/r4.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(31, 45.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext4.2", log = F)

writeRaster(r4.3_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext4.3/r4.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(46, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext4.3", log = F)

writeRaster(r5.1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext5.1/r5.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(61, 90.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext5.1", log = F)

writeRaster(r5.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext5.2/r5.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(91, 105.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext5.2", log = F)

writeRaster(r5.3_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext5.3/r5.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(106, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext5.3", log = F)

writeRaster(r6.1.1.1_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.1.1/r6.1.1.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(121, 128.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.1.1", log = F)

writeRaster(r6.1.1.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.1.2/r6.1.1.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(129, 135.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.1.2", log = F)

writeRaster(r6.1.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.2/r6.1.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(136, 150.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext6.1.2", log = F)

writeRaster(r6.2_sq, "/Volumes/Personal/lab/raster/squamata_svl/ext6.2/r6.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_svl", 
        extent(151, 180, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_svl/ext6.2", log = F)


writeRaster(r_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/r.tif",
            format = "GTiff")

writeRaster(r1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext1/r1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-180, -120, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext1", log = F)

writeRaster(r2.1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.1/r2.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-119.9999, -100, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.1", log = F)

writeRaster(r2.2.1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.2.1/r2.2.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-99.9999, -85, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.2.1", log = F)

writeRaster(r2.2.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.2.2/r2.2.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-84.9999, -75, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.2.2", log = F)

writeRaster(r2.3_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.3/r2.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-74.9999, -40, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext2.3", log = F)

writeRaster(r3_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext3/r3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(-39.9999, 0.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext3", log = F)

writeRaster(r4.1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.1/r4.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(1, 30.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.1", log = F)

writeRaster(r4.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.2/r4.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(31, 45.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.2", log = F)

writeRaster(r4.3_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.3/r4.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(46, 60.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext4.3", log = F)

writeRaster(r5.1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.1/r5.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(61, 90.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.1", log = F)

writeRaster(r5.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.2/r5.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(91, 105.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.2", log = F)

writeRaster(r5.3_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.3/r5.3.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(106, 120.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext5.3", log = F)

writeRaster(r6.1.1.1_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.1.1/r6.1.1.1.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(121, 128.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.1.1", log = F)

writeRaster(r6.1.1.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.1.2/r6.1.1.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(129, 135.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.1.2", log = F)

writeRaster(r6.1.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.2/r6.1.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(136, 150.9999, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.1.2", log = F)

writeRaster(r6.2_sq, "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.2/r6.2.tif",
            format = "GTiff")
crop_rg("/Volumes/Personal/lab/raster/squamata_rates_svl", 
        extent(151, 180, -90, 90),
        "/Volumes/Personal/lab/raster/squamata_rates_svl/ext6.2", log = F)

get_values <- function(full) {

  data(wrld_simpl, package = "maptools")

  rem1 <- extract(full, wrld_simpl, cellnumbers = T, weights = T, 
                  small = T)
  rem2 <- do.call(rbind.data.frame, rem1)[, 1]
  values(full)[-rem2] <- NA

  brks <- quantile(values(full)[order(values(full))], probs = seq(0, 1, 0.02),
                   na.rm = T)
  brks <- unique(brks)
  full_p <- rasterToPoints(full, spatial = TRUE)
  full_df  <- data.frame(full_p)
  colnames(full_df)[1] <- "index_1"
  full_df <- full_df %>% mutate(index_2 = cut(index_1, breaks = brks))

  return(full_df)
}

raster_svl_am <- get_values(raster("data/amphibia/svl.tif"))
raster_rates2_am <- get_values(raster("data/amphibia/rates_svl.tif"))

raster_svl_sq <- get_values(raster("data/reptilia/svl.tif"))
raster_rates2_sq <- get_values(raster("data/reptilia/rates_svl.tif"))

r <- raster(ncols = 2160, nrows = 900, ymn = -60)

colors_svl <- mako(50)[50:1]
colors_rates <- rocket(50)[50:1]

ras_wrld <- rasterToPoints(rasterize(wrld_simpl, r), spatial = TRUE)
ras_wrld_df  <- data.frame(ras_wrld)

plot1_svl_am <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_svl_am, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_svl) +
  theme_void() 

plot2_svl_am <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates2_am, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

plot1_svl_sq <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_svl_sq, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_svl) +
  theme_void() 

plot2_svl_sq <- ggplot() +
  geom_raster(data = ras_wrld_df, aes(x=x, y=y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_rates2_sq, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_rates) +
  theme_void() 

inches <- 4.5
res <- 600

tiff("figures/FigureS8_am1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_svl_am)  
dev.off() 

tiff("figures/FigureS8_am2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_svl_am)  
dev.off() 

tiff("figures/FigureS8_sq1.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot1_svl_sq)  
dev.off() 

tiff("figures/FigureS8_sq2.tiff", width = inches*res, height = inches*res/2, 
     units = "px")
print(plot2_svl_sq)  
dev.off() 

quantile(round(unique(quantile(values(raster("data/amphibia/svl.tif"))[order(values(raster("data/amphibia/svl.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/amphibia/rates_svl.tif"))[order(values(raster("data/amphibia/rates_svl.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/reptilia/svl.tif"))[order(values(raster("data/reptilia/svl.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

quantile(round(unique(quantile(values(raster("data/reptilia/rates_svl.tif"))[order(values(raster("data/reptilia/rates_svl.tif")))], probs = seq(0, 1, 0.02), na.rm = T)), 3))

