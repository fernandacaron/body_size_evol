rm(list  =  ls())

setwd("~/Documents/lab/body_size_evol")

library(viridis)
library(phytools)
library(geiger)
library(moments)

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_RMA_17jan24.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

colorRampAlpha  <-  function(..., n, alpha) {
  colors  <-  colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep = "")
}

col_al_am <- colorRampAlpha(c("plum", "orchid4", "black"), n = 3, alpha = 0.6)
col_am <- colorRampAlpha(c("plum", "orchid4", "black"), n = 3, alpha = 1)

col_al_sq <- colorRampAlpha(c("lightsteelblue1", "slateblue2", "black"), n = 6,
                            alpha = 0.6)
col_sq <- colorRampAlpha(c("lightsteelblue1", "slateblue2", "black"), n = 6,
                         alpha = 1)

col_al_av <- colorRampAlpha(c("darkolivegreen1", "black"), n = 7, alpha = 0.6)
col_av <- colorRampAlpha(c("darkolivegreen1", "black"), n = 7, alpha = 1)

col_al_ma <- colorRampAlpha(c("lightgoldenrod1", "orange", "black"), n = 7,
                            alpha = 0.6)
col_ma <- colorRampAlpha(c("lightgoldenrod1", "orange", "black"), n = 7,
                         alpha = 1)

cols_alp <- c(rgb(139/255,71/255,137/255,0.4), rgb(97/255,82/255,190/255,0.4),
           rgb(161/255,204/255,89/255,0.5), rgb(255/255,207/255,83/255,0.5))

cols <- c(rgb(139/255,71/255,137/255), rgb(97/255,82/255,190/255),
        rgb(161/255,204/255,89/255), rgb(255/255,207/255,83/255))

#Figure 1 - Histogramas log body size usando as m?dias das fontes 
pdf("EvolBiol_20oct23/rev1/Figure1.pdf", width = 10, height = 9)

#Amphibia - Gimnophiona n?o, s? 8 spp

par(fig = c(0, 0.5, 0.66, 1), mar = c(2, 4, 2, 2))

hist(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Anura"]), ylim = c(0, 1.2), 
     main = "", xlim = c(-3, 15),  xlab = "", breaks = 20, col = col_al_am[1], 
     freq = F, border = col_am[1])
hist(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Caudata"]), add = T, 
     col = col_al_am[2], breaks = 20, freq = F, border = col_am[2])
hist(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Gymnophiona"]), add = T, 
     col = col_al_am[3], breaks = 10, freq = F, border = col_am[3])
title("A", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_am, 
       legend = c("Anura", "Caudata", "Gymnophiona"))

skewness(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Anura"]), na.rm = T)
skewness(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Caudata"]), na.rm = T)
skewness(log(dat_am$Body_mass_g_RMA[dat_am$Order == "Gymnophiona"]), na.rm = T)

#Squamata

par(fig = c(0.5, 1, 0.66, 1), mar = c(2, 3, 2, 2), new = T)

hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Anguimorpha"]),
     ylim = c(0, 0.5), xlim = c(-3, 15), main = "", xlab = "", ylab = "", 
     breaks = 20, col = col_al_sq[1], border = col_sq[1], freq = F)
hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Gekkota"]), add = T, 
     breaks = 20, col = col_al_sq[2], freq = F, border = col_sq[2])
hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Iguania"]), add = T, 
     breaks = 20, col = col_al_sq[3], freq = F, border = col_sq[3])
hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Lacertoidea"]), add = T, 
     breaks = 20, col = col_al_sq[4], freq = F, border = col_sq[4])
hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Scincoidea"]), add = T,
     breaks = 20, col = col_al_sq[5], freq = F, border = col_sq[5])
hist(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Serpentes"]), add = T, 
     breaks = 20, col = col_al_sq[6], freq = F, border = col_sq[6])
title("B", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_sq, 
       legend = c("Anguimorpha",  "Gekkota", "Iguania", "Lacertoidea", 
                  "Scincoidea", "Serpentes"))

skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Anguimorpha"]), na.rm = T)
skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Gekkota"]), na.rm = T)
skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Iguania"]), na.rm = T)
skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Lacertoidea"]), na.rm = T)
skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Scincoidea"]), na.rm = T)
skewness(log(dat_sq$Body_mass_g_mean[dat_sq$Taxon == "Serpentes"]), na.rm = T)

#Aves

par(fig = c(0, 0.5, 0.33, 0.66), mar = c(2 ,4, 2, 0), new = T)

hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Apodiformes"]), 
     ylim = c(0, 1), xlim = c(0, 10), main = "", xlab = "", breaks = 30, 
     col = col_al_av[1], freq = F, border = col_av[1])
hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Charadriiformes"]), add = T, 
     breaks = 20, col = col_al_av[2], freq = F, border = col_av[2])
hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Columbiformes"]), add = T, 
     breaks = 20, col = col_al_av[3], freq = F, border = col_av[3])
hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Passeriformes"]), add = T, 
     breaks = 30, col = col_al_av[4], freq = F, border = col_av[4])
hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Piciformes"]), add = T, 
     breaks = 20, col = col_al_av[5], freq = F, border = col_av[5])
hist(log(dat_av$Body_mass_g_mean[dat_av$Order == "Psittaciformes"]), add = T,
     breaks = 20, col = col_al_av[6], freq = F, border = col_av[6])
hist(log(dat_av$Body_mass_g_mean[!dat_av$Order %in% c("Apodiformes",
                                                      "Charadriiformes",
                                                      "Columbiformes",
                                                      "Passeriformes",
                                                      "Piciformes",
                                                      "Psittaciformes")]), add = T,
     breaks = 40, col = col_al_av[7], freq = F, border = col_av[7])
title("C", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_av, 
       legend  = c("Apodiformes", "Charadriiformes", "Columbiformes",
                   "Passeriformes", "Piciformes", "Psittaciformes", "Other"))

skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Apodiformes"]), na.rm = T)
skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Charadriiformes"]), na.rm = T)
skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Columbiformes"]), na.rm = T)
skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Passeriformes"]), na.rm = T)
skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Piciformes"]), na.rm = T)
skewness(log(dat_av$Body_mass_g_mean[dat_av$Order == "Psittaciformes"]), na.rm = T)


#Mammalia

par(fig = c(0.5, 1, 0.33, 0.66), mar = c(2, 3, 2, 2), new = T)

hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Carnivora"]), ylim = c(0, 1),	
     xlim = c(0, 21), main = "", xlab = "", ylab = "", breaks = 20, 
     col = col_al_ma[1], freq = F, border = col_ma[1])
hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Cetartiodactyla"]), add = T, 
     breaks = 20, col = col_al_ma[2], freq = F, border = col_ma[2])
hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Chiroptera"]), add = T, 
     breaks = 20, col = col_al_ma[3], freq = F, border = col_ma[3])
hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Eulipotyphla"]), add = T, 
     breaks = 20, col = col_al_ma[4], freq = F, border = col_ma[4])
hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Primates"]), add = T, 
     breaks = 20, col = col_al_ma[5], freq = F, border = col_ma[5])
hist(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Rodentia"]), add = T, 
     breaks = 20, col = col_al_ma[6], freq = F, border = col_ma[6])
hist(log(dat_ma$Body_mass_g_mean[!dat_ma$Order %in% c("Carnivora",
                                                      "Cetartiodactyla",
                                                      "Chiroptera",
                                                      "Eulipotyphla",
                                                      "Primates",
                                                      "Rodentia")]), add = T,
     breaks = 40, col = col_al_ma[7], freq = F, border = col_ma[7])
title("D", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_ma, 
       legend = c("Carnivora", "Cetartiodactyla", "Chiroptera", "Eulipotyphla",
                  "Primates", "Rodentia", "Other"))

skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Carnivora"]), na.rm = T)
skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Cetartiodactyla"]), na.rm = T)
skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Chiroptera"]), na.rm = T)
skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Eulipotyphla"]), na.rm = T)
skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Primates"]), na.rm = T)
skewness(log(dat_ma$Body_mass_g_mean[dat_ma$Order == "Rodentia"]), na.rm = T)

#Todos juntos

par(fig = c(0, 1, 0, 0.33), mar = c(5, 4, 2, 2), new = T)

hist(log(dat_am$Body_mass_g_RMA), ylim = c(0,0.9), xlim = c(-5,20), freq = F, 
     breaks = 60, main = "", xlab = "log Body Mass (g)", col = cols_alp[1], 
     border = cols[1])
abline(v = median(log(dat_am$Body_mass_g_RMA), na.rm = T), col = cols[1], lwd = 2)

hist(log(dat_sq$Body_mass_g_mean), add = T, breaks = 60, freq = F, 
     col = cols_alp[2], border = cols[2])
abline(v = median(log(dat_sq$Body_mass_g_mean), na.rm = T), col = cols[2], 
       lwd = 2)

hist(log(dat_av$Body_mass_g_mean), add = T, breaks = 60, freq = F, 
     col = cols_alp[3],border = cols[3])
abline(v = median(log(dat_av$Body_mass_g_mean), na.rm = T), col = cols[3], 
       lwd = 2)

hist(log(dat_ma$Body_mass_g_mean), add = T, breaks = 70, freq = F, 
     col = cols_alp[4],
     border = cols[4])
abline(v = median(log(dat_ma$Body_mass_g_mean), na.rm = T), col = cols[4], 
       lwd = 2)

title("E", adj = 0)
legend("topright", pch = 15, bty = 'n', col = cols_alp, 
       legend = c("Amphibia", "Squamata", "Aves", "Mammalia"))

dev.off()

#################

#Figure S1 - Median
pdf("EvolBiol_20oct23/rev1/FigureS1.pdf", width  =  14)

#Squamata

par(fig = c(0, 0.33, 0.5, 1), mar = c(2, 3, 2, 1))

hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Anguimorpha"]), 
     ylim = c(0, 0.9), xlim = c(-4, 15), main = "", xlab = "", ylab = "", 
     breaks = 20, col = col_al_sq[1], border = col_sq[1], freq = F)
hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Gekkota"]), add = T, 
     breaks = 20, col = col_al_sq[2], freq = F, border = col_sq[2])
hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Iguania"]), add = T, 
     breaks = 20, col = col_al_sq[3], freq = F, border = col_sq[3])
hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Lacertoidea"]), add = T, 
     breaks = 20, col = col_al_sq[4], freq = F, border = col_sq[4])
hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Scincoidea"]), add = T, 
     breaks = 20, col = col_al_sq[5], freq = F, border = col_sq[5])
hist(log(dat_sq$Body_mass_g_median[dat_sq$Taxon == "Serpentes"]), add = T, 
     breaks = 20, col = col_al_sq[6], freq = F, border = col_sq[6])
title("A", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_sq, 
       legend = c("Anguimorpha",  "Gekkota", "Iguania", "Lacertoidea", 
                  "Scincoidea", "Serpentes"))
				
#Aves

par(fig = c(0.33, 0.66, 0.5, 1), mar = c(2, 3, 2, 1), new = T)
 
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Apodiformes"]), 
     ylim = c(0, 0.9), xlim = c(0, 11), main = "", xlab = "", breaks = 30, 
     col = col_al_av[1], freq = F, border = col_av[1])
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Charadriiformes"]), add = T,
     breaks = 20, col = col_al_av[2], freq = F, border = col_av[2])
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Columbiformes"]), add = T, 
     breaks = 20, col = col_al_av[3], freq = F, border = col_av[3])
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Passeriformes"]), add = T, 
     breaks = 30, col = col_al_av[4], freq = F, border = col_av[4])
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Piciformes"]), add = T, 
     breaks = 20, col = col_al_av[5], freq = F, border = col_av[5])
hist(log(dat_av$Body_mass_g_median[dat_av$Order == "Psittaciformes"]), add = T,
     breaks = 20, col = col_al_av[6], freq = F, border = col_av[6])
hist(log(dat_av$Body_mass_g_median[!dat_av$Order %in% c("Apodiformes",
                                                      "Charadriiformes",
                                                      "Columbiformes",
                                                      "Passeriformes",
                                                      "Piciformes",
                                                      "Psittaciformes")]), add = T,
     breaks = 40, col = col_al_av[7], freq = F, border = col_av[7])
title("B", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_av, 
       legend  = c("Apodiformes", "Charadriiformes", "Columbiformes",
                   "Passeriformes", "Piciformes", "Psittaciformes", "Other"))

#Mammalia

par(fig = c(0.66 , 1, 0.5 ,1), mar = c(2, 3, 2, 1), new = T)

hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Carnivora"]), 
     ylim = c(0, 0.9), xlim = c(0, 21), main = "", xlab = "", ylab = "", 
     breaks = 20, col = col_al_ma[1], freq = F, border = col_ma[1])
hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Cetartiodactyla"]), add = T,
     breaks = 20, col = col_al_ma[2], freq = F, border = col_ma[2])
hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Chiroptera"]), add = T, 
     breaks = 20, col = col_al_ma[3], freq = F, border = col_ma[3])
hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Eulipotyphla"]), add = T, 
     breaks = 20, col = col_al_ma[4], freq = F, border = col_ma[4])
hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Primates"]), add = T, 
     breaks = 20, col = col_al_ma[5], freq = F, border = col_ma[5])
hist(log(dat_ma$Body_mass_g_median[dat_ma$Order == "Primates"]), add = T, 
     breaks = 20, col = col_al_ma[6], freq = F, border = col_ma[6])
hist(log(dat_ma$Body_mass_g_median[!dat_ma$Order %in% c("Carnivora",
                                                      "Cetartiodactyla",
                                                      "Chiroptera",
                                                      "Eulipotyphla",
                                                      "Primates",
                                                      "Rodentia")]), add = T,
     breaks = 40, col = col_al_ma[7], freq = F, border = col_ma[7])
title("C", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_ma, 
       legend = c("Carnivora", "Cetartiodactyla", "Chiroptera", "Eulipotyphla",
                  "Primates", "Rodentia", "Other"))


#Todos juntos

par(fig = c(0, 1, 0, 0.5), mar = c(5, 4, 2, 2), new = T)

hist(log(dat_sq$Body_mass_g_median), ylim = c(0,0.32), xlim = c(-5, 20), 
     freq = F, breaks = 60, main = "", xlab = "log Body Mass (g)", 
     col = cols_alp[2], border = cols[2])
abline(v = median(log(dat_sq$Body_mass_g_median), na.rm = T), col = cols[1], 
       lwd = 2)

hist(log(dat_av$Body_mass_g_median), add = T, breaks = 60, freq = F, 
     col = cols_alp[3], border = cols[3])
abline(v = median(log(dat_av$Body_mass_g_median), na.rm = T), col = cols[3], 
       lwd = 2)

hist(log(dat_ma$Body_mass_g_median), add = T, breaks = 70, freq = F, 
     col = cols_alp[4], border = cols[4])
abline(v = median(log(dat_ma$Body_mass_g_median), na.rm = T), col = cols[4], 
       lwd = 2)

title("D", adj = 0)
legend("topright", pch = 15, bty = 'n', col = cols_alp[2:4], 
       legend = c("Squamata", "Aves", "Mammalia"))

dev.off()

##################

#Figure S2 - Maximum
pdf("EvolBiol_20oct23/rev1/FigureS2.pdf", width  =  14)

#Squamata

par(fig = c(0, 0.33, 0.5, 1), mar = c(2, 3, 2, 1))

hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Anguimorpha"]), ylim = c(0, 0.9),
     xlim = c(-4, 15), main = "", xlab = "", ylab = "", breaks = 20, 
     col = col_al_sq[1], border = col_sq[1], freq = F)
hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Gekkota"]), add = T, breaks = 20,
     col = col_al_sq[2], freq = F, border = col_sq[2])
hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Iguania"]), add = T, breaks = 20,
     col = col_al_sq[3], freq = F, border = col_sq[3])
hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Lacertoidea"]), add = T, 
     breaks = 20, col = col_al_sq[4], freq = F, border = col_sq[4])
hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Scincoidea"]), add = T, 
     breaks = 20, col = col_al_sq[5], freq = F, border = col_sq[5])
hist(log(dat_sq$Body_mass_g_max[dat_sq$Taxon == "Serpentes"]), add = T, 
     breaks = 40, col = col_al_sq[6], freq = F, border = col_sq[6])
title("A", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_sq, 
       legend = c("Anguimorpha",  "Gekkota", "Iguania", "Lacertoidea", 
                  "Scincoidea", "Serpentes"))

#Aves

par(fig = c(0.33, 0.66, 0.5, 1), mar = c(2, 3, 2, 1), new = T)

hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Apodiformes"]), 
     ylim = c(0,0.9), xlim = c(0, 11), main = "", xlab = "", breaks = 30, 
     col = col_al_av[1], freq = F, border = col_av[1])
hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Charadriiformes"]), add = T, 
     breaks = 20, col = col_al_av[2], freq = F, border = col_av[2])
hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Columbiformes"]), add = T, 
     breaks = 20, col = col_al_av[3], freq = F, border = col_av[3])
hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Passeriformes"]), add = T, 
     breaks = 30, col = col_al_av[4], freq = F, border = col_av[4])
hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Piciformes"]), add = T, 
     breaks = 20, col = col_al_av[5], freq = F, border = col_av[5])
hist(log(dat_av$Body_mass_g_max[dat_av$Order == "Psittaciformes"]), add = T,
     breaks = 20, col = col_al_av[6], freq = F, border = col_av[6])
hist(log(dat_av$Body_mass_g_max[!dat_av$Order %in% c("Apodiformes",
                                                        "Charadriiformes",
                                                        "Columbiformes",
                                                        "Passeriformes",
                                                        "Piciformes",
                                                        "Psittaciformes")]), add = T,
     breaks = 40, col = col_al_av[7], freq = F, border = col_av[7])
title("B", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_av, 
       legend  = c("Apodiformes", "Charadriiformes", "Columbiformes",
                   "Passeriformes", "Piciformes", "Psittaciformes", "Other"))

#Mammalia

par(fig = c(0.66, 1, 0.5, 1), mar = c(2, 3, 2, 1), new = T)

hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Carnivora"]), ylim = c(0, 0.9),
     xlim = c(0, 21), main = "", xlab = "", ylab = "", breaks = 20, 
     col = col_al_ma[1], freq = F, border = col_ma[1])
hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Cetartiodactyla"]), add = T, 
     breaks = 20, col = col_al_ma[2], freq = F, border = col_ma[2])
hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Chiroptera"]), add = T, 
     breaks = 20, col = col_al_ma[3], freq = F, border = col_ma[3])
hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Eulipotyphla"]), add = T, 
     breaks = 20, col = col_al_ma[4], freq = F, border = col_ma[4])
hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Primates"]), add = T, 
     breaks = 20, col = col_al_ma[5], freq = F, border = col_ma[5])
hist(log(dat_ma$Body_mass_g_max[dat_ma$Order == "Rodentia"]), add = T, 
     breaks = 20, col = col_al_ma[6], freq = F, border = col_ma[6])
hist(log(dat_ma$Body_mass_g_max[!dat_ma$Order %in% c("Carnivora",
                                                        "Cetartiodactyla",
                                                        "Chiroptera",
                                                        "Eulipotyphla",
                                                        "Primates",
                                                        "Rodentia")]), add = T,
     breaks = 40, col = col_al_ma[7], freq = F, border = col_ma[7])
title("C", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_al_ma, 
       legend = c("Carnivora", "Cetartiodactyla", "Chiroptera", "Eulipotyphla",
                  "Primates", "Rodentia", "Other"))


#Todos juntos

par(fig = c(0, 1, 0, 0.5), mar = c(5, 4, 2, 2), new = T)

hist(log(dat_sq$Body_mass_g_max), ylim = c(0, 0.32), xlim = c(-5, 21), freq = F,
     main = "", xlab = "log Body Mass (g)", breaks = 60, col = cols_alp[2], 
     border = cols[2])
abline(v = median(log(dat_sq$Body_mass_g_max), na.rm  =  T), col = cols[2], 
       lwd = 2)

hist(log(dat_av$Body_mass_g_max), add = T, breaks = 60, freq = F, 
     col = cols_alp[3], border = cols[3])
abline(v = median(log(dat_av$Body_mass_g_max), na.rm  =  T), col = cols[3], 
       lwd = 2)

hist(log(dat_ma$Body_mass_g_max), add = T, breaks = 70, freq = F, 
     col = cols_alp[4], border = cols[4])
abline(v = median(log(dat_ma$Body_mass_g_max), na.rm  =  T), col = cols[4], 
       lwd = 2)

title("D", adj = 0)
legend("topright", pch = 15, bty = 'n', col = cols_alp[2:4], 
       legend = c("Squamata", "Aves", "Mammalia"))

dev.off()

###############

#Figure 2 - mapping body size in the phylogenies
tr_am <- read.nexus("~/Documents/lab/data/trees/amphibia_VertLife_27JUL20.nex")
tr_sq <- read.nexus("~/Documents/lab/data/trees/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("~/Documents/lab/data/trees/mammalia_node_dated_VertLife_27JUL20.nex")

log_mass_am <- log(dat_am$Body_mass_g_RMA)
names(log_mass_am) <- dat_am$Scientific_name

log_mass_sq <- log(dat_sq$Body_mass_g_mean)
names(log_mass_sq) <- dat_sq$Scientific_name

log_mass_av <- log(dat_av$Body_mass_g_mean)
names(log_mass_av) <- dat_av$Scientific_name

log_mass_ma <- log(dat_ma$Body_mass_g_mean)
names(log_mass_ma) <- dat_ma$Scientific_name

log_mass_am <- log_mass_am[complete.cases(log_mass_am)]
log_mass_sq <- log_mass_sq[complete.cases(log_mass_sq)]
log_mass_av <- log_mass_av[complete.cases(log_mass_av)]
log_mass_ma <- log_mass_ma[complete.cases(log_mass_ma)]

tr_map_am <- treedata(tr_am[[1]], log_mass_am)$phy
log_mass_am <- log_mass_am[tr_map_am$tip.label]

tr_map_sq <- treedata(tr_sq[[1]], log_mass_sq)$phy
log_mass_sq <- log_mass_sq[tr_map_sq$tip.label]

tr_map_av <- treedata(tr_av[[1]], log_mass_av)$phy
log_mass_av <- log_mass_av[tr_map_av$tip.label]

tr_map_ma <- treedata(tr_ma[[1]], log_mass_ma)$phy
log_mass_ma <- log_mass_ma[tr_map_ma$tip.label]

cols_pal <- turbo(20)

map_am <- contMap(tr_map_am, log_mass_am, plot = F, lims = c(-2.5, 19))
plot_am <- setMap(map_am, cols_pal)

map_sq <- contMap(tr_map_sq, log_mass_sq, plot = F, lims = c(-2.5, 19))
plot_sq <- setMap(map_sq, cols_pal)

map_av <- contMap(tr_map_av, log_mass_av, plot = F, lims = c(-2.5, 19))
plot_av <- setMap(map_av, cols_pal)

map_ma <- contMap(tr_map_ma, log_mass_ma, plot = F, lims = c(-2.5, 19))
plot_ma <- setMap(map_ma, cols_pal)

pdf("EvolBiol_20oct23/rev1/Figure2_median.pdf", width = 8)

layout(matrix(c(1,2,3,4,5,5), nrow = 2), widths  =  c(1, 1, 0.4))

par(mar = c(3,4,4,4))

plot(plot_am, ftype = "off", type = "fan", lwd = 1, outline = F, fsize = 1, legend = F,
     mar = c(2,1,1,1))
title("A", adj = 0, line = -1)
mtext("Amphibia", cex = 0.6, side = 1)

plot(plot_av, ftype = "off", type = "fan", lwd = 1, outline = F, fsize = 1, legend = F, mar = c(2,1,1,1))
title("C", adj = 0, line = -1)
mtext("Aves", cex = 0.6, side = 1)

plot(plot_sq, ftype = "off", type = "fan", lwd = 1, outline = F, fsize = 1, legend = F, mar = c(2,1,1,1))
title("B", adj = 0, line = -1)
mtext("Squamata", cex = 0.6, side = 1)

plot(plot_ma, ftype = "off", type = "fan", lwd = 1, outline = F, fsize = 1, legend = F,  mar = c(2,1,1,1))
title("D", adj = 0, line = -1)
mtext("Mammalia", cex = 0.6, side = 1)

par(mar = c(4, 4, 3, 4))

color.bar  <-  function(lut, min, max  =  -min, nticks  =  11, 
                      ticks  =  seq(min, max, len  =  nticks), title  =  "") {
  scale  =  (length(lut) - 1)/(max - min)
  
  plot(c(0, 10), c(min, max), type  =  "n", bty  =  "n", xaxt  =  "n", 
       xlab  =  "", yaxt  =  "n", ylab  =  "", main  =  title)
  axis(2, ticks, las  =  1)
  for (i in 1:(length(lut) - 1)) {
    y  =  (i - 1)/scale + min
    rect(0, y, 10, y + 1/scale, col  =  lut[i], border  =  NA)
  }
}

color.bar(turbo(200), min  =  -2.5, max  =  19, title  =  "")
mtext("log Body Size (g)", cex = 0.6)

dev.off()
