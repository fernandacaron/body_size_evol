rm(list = ls())

setwd("~/Documents/lab/body_size_evol")

library(phytools)
library(geiger)
library(RRphylo)
library(raster)
library(stringi)
library(fasterize)
library(sf)

tr_am <- read.nexus("~/Documents/lab/data/trees/amphibia_VertLife_27JUL20.nex")
tr_sq <- read.nexus("~/Documents/lab/data/trees/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("~/Documents/lab/data/trees/mammalia_node_dated_VertLife_27JUL20.nex")

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

mass_am <- log(dat_am$Body_mass_g_1[complete.cases(dat_am$Body_mass_g_1)])
names(mass_am) <- dat_am$Scientific_name[complete.cases(dat_am$Body_mass_g_1)]

mass_sq <- log(dat_sq$Body_mass_g_mean[complete.cases(dat_sq$Body_mass_g_mean)])
names(mass_sq) <- dat_sq$Scientific_name[complete.cases(dat_sq$Body_mass_g_mean)]

mass_av <- log(dat_av$Body_mass_g_mean[complete.cases(dat_av$Body_mass_g_mean)])
names(mass_av) <- dat_av$Scientific_name[complete.cases(dat_av$Body_mass_g_mean)]

mass_ma <- log(dat_ma$Body_mass_g_mean[complete.cases(dat_ma$Body_mass_g_mean)])
names(mass_ma) <- dat_ma$Scientific_name[complete.cases(dat_ma$Body_mass_g_mean)]

RR_am <- list()
for (i in 1:100) {
  print(i)
  
  pruned_tr_am <- treedata(tr_am[[i]], mass_am, warnings = F)$phy
  
  mass_am <- mass_am[names(mass_am) %in% pruned_tr_am$tip.label]
  
  mass_am <- mass_am[match(pruned_tr_am$tip.label, names(mass_am))]
  
  RR_am[[i]] <- RRphylo(tree = pruned_tr_am, y = mass_am)
}
#save(RR_am, file = "data/amphibia/RR_am.RData")
#load(file = "data/amphibia/RR_am.RData")

RR_sq <- list()
for (i in 1:100) {
  print(i)
  
  pruned_tr_sq <- treedata(tr_sq[[i]], mass_sq, warnings = F)$phy
  
  mass_sq <- mass_sq[names(mass_sq) %in% pruned_tr_sq$tip.label]
  
  mass_sq <- mass_sq[match(pruned_tr_sq$tip.label, names(mass_sq))]
  
  RR_sq[[i]] <- RRphylo(tree = pruned_tr_sq, y = mass_sq)
}
#save(RR_sq, file = "data/squamata/RR_sq.RData")
#load(file = "data/squamata/RR_sq.RData")

RR_av <- list()
for (i in 1:100) {
  print(i)
  
  pruned_tr_av <- treedata(tr_av[[i]], mass_av, warnings = F)$phy
  
  mass_av <- mass_av[names(mass_av) %in% pruned_tr_av$tip.label]
  
  mass_av <- mass_av[match(pruned_tr_av$tip.label, names(mass_av))]
  
  RR_av[[i]] <- RRphylo(tree = pruned_tr_av, y = mass_av)
}
#save(RR_av, file = "data/aves/RR_av.RData")
#load(file = "data/aves/RR_av.RData")

RR_ma <- list()
for (i in 1:100) {
  print(i)
  
  pruned_tr_ma <- treedata(tr_ma[[i]], mass_ma, warnings = F)$phy
  
  mass_ma <- mass_ma[names(mass_ma) %in% pruned_tr_ma$tip.label]
  
  mass_ma <- mass_ma[match(pruned_tr_ma$tip.label, names(mass_ma))]
  
  RR_ma[[i]] <- RRphylo(tree = pruned_tr_ma, y = mass_ma)
}
#save(RR_ma, file = "data/mammalia/RR_ma.RData")
#load(file = "data/mammalia/RR_ma.RData")


RR_am_res <- matrix(nrow = length(mass_am), ncol = 10)
rownames(RR_am_res) <- pruned_tr_am$tip.label[order(pruned_tr_am$tip.label)]

for (i in 1:length(mass_am)) {
  spp <- rownames(RR_am_res)[i]
  for (j in 1:10) {
    RR_am_res[i, j] <- RR_am[[j]]$rates[rownames(RR_am[[j]]$rates) == spp]
  }
}

RR_am_mean <- rowMeans(RR_am_res)
names(RR_am_mean) <- rownames(RR_am_res)

map_am <- shapefile("/Users/fercaron/Documents/lab/arrested_div/data/AMPHIBIANS_IUCN/AMPHIBIANS.shp")
map_am$binomial <- stri_replace_all_fixed(map_am$binomial, " ", "_")
map_am2 <- raster::subset(map_am, map_am$binomial %in% rownames(RR_am_res))

r <- raster(ncols = 2160, nrows = 900, ymn = -60)
raster_stack <- r

RR_am_mean <- RR_am_mean[names(RR_am_mean) %in% map_am2$binomial]

for (i in 1:length(RR_am_mean)) {
  print(i) 
  
  s <- as.character(names(RR_am_mean)[i]) 
  map_i <- subset(map_am2, map_am2$binomial == s)
  
  raster_i <- rasterize(map_i, r)
  
  rastercells <- which(getValues(!is.na(raster_i))) 
  raster_i[rastercells] <- RR_am_mean[i]
  
  if(i == 1) {
    raster_stack <- raster_i
  } else {
    raster_stack <- addLayer(raster_stack, raster_i)
  }
}

full_am <- stackApply(raster_stack, indices = rep(1, length(RR_am_mean)), 
                          fun = median, na.rm = T) 

for (i in 1:length(sdi_mal)) {
  print(i) 
  
  s <- as.character(names(sdi_mal)[i]) 
  map_i <- subset(birds2, birds2$sci_name == s)
  
  for (j in 1:length(map_i$Shape)) { 
    try(map_i$Shape[[j]] <- st_cast(map_i$Shape[[j]], 'MULTIPOLYGON'))
  }
  try(map_i$Shape <- st_cast(map_i$Shape, 'MULTIPOLYGON'))
  
  raster_i <- fasterize(st_as_sf(map_i$Shape), r)
  
  rastercells <- which(getValues(!is.na(raster_i))) 
  raster_i[rastercells] <- sdi_mal[i]
  
  if(i == 1) {
    raster_stack_mal <- raster_i
  } else {
    raster_stack_mal <- addLayer(raster_stack_mal, raster_i)
  }
}

av_full_mal <- stackApply(raster_stack_mal, indices = rep(1, length(sdi_mal)), 
                          fun = median, na.rm = T) 

data(wrld_simpl)
rem_fem1 <- extract(av_full_fem, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_fem2 <- do.call(rbind.data.frame, rem_fem1)[, 1]
values(av_full_fem)[-rem_fem2] <- NA

rem_mal1 <- extract(av_full_mal, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_mal2 <- do.call(rbind.data.frame, rem_mal1)[, 1]
values(av_full_mal)[-rem_mal2] <- NA

brks_fem <- quantile(values(av_full_fem)[order(values(av_full_fem))], 
                     probs = seq(0, 1, 0.02), na.rm = T)
av_full_fem_p <- rasterToPoints(av_full_fem, spatial = TRUE)
av_full_fem_df  <- data.frame(av_full_fem_p)
av_full_fem_df <- av_full_fem_df %>% mutate(index_2 = cut(index_1, 
                                                          breaks = brks_fem))

brks_mal <- quantile(values(av_full_mal)[order(values(av_full_mal))], 
                     probs = seq(0, 1, 0.02), na.rm = T)
av_full_mal_p <- rasterToPoints(av_full_mal, spatial = TRUE)
av_full_mal_df  <- data.frame(av_full_mal_p)
av_full_mal_df <- av_full_mal_df %>% mutate(index_2 = cut(index_1, 
                                                          breaks = brks_mal))

colors_fem <- met.brewer(name = "Hiroshige", n = 50, direction = -1)[1:50]
colors_mal <- met.brewer(name = "Hiroshige", n = 50, direction = 1)[1:50]

inches <- 4.5
res <- 600

ggplot_fem <- ggplot() +
  geom_raster(data = av_full_fem_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_fem) +
  theme_void() 

plot_fem <- "Figure3_female.tiff"
tiff(plot_fem, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_fem)  
dev.off()  

ggplot_mal <- ggplot() +
  geom_raster(data = av_full_mal_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mal) +
  theme_void() 

plot_mal <- "Figure3_male.tiff"
tiff(plot_mal, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_mal)  
dev.off() 

plot_sca<- "Figure3_ScaleBar.tiff"
tiff(plot_sca, width = 1*res, h = 0.1*res, units = "px")
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
      col = colors_fem, breaks = breaks)
axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                  cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                  labels = rep("", 5), tck = 0.2)
do.call("axis", axis.args)
box()
dev.off()

as.numeric(quantile(values(av_full_fem)[order(values(av_full_fem))], 
                    probs = seq(0, 1, 0.01), na.rm = T)[seq(1, 101, 25)])

as.numeric(quantile(values(av_full_mal)[order(values(av_full_mal))], 
                    probs = seq(0, 1, 0.01), na.rm = T)[seq(1, 101, 25)])
