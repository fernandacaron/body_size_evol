rm(list = ls())

setwd("~/Documents/lab/body_size_evol")

library(geiger)
library(GET)

## Carregando os dados
tr_am <- read.nexus("~/Documents/lab/data/trees/amphibia_VertLife_27JUL20.nex")[1:100]
tr_sq <- read.nexus("~/Documents/lab/data/trees/squamata_VertLife_27JUL20.nex")[1:100]
tr_av <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")[1:100]
tr_ma <- read.nexus("~/Documents/lab/data/trees/mammalia_node_dated_VertLife_27JUL20.nex")[1:100]

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_RMA_17jan24.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

mass_am <- log(dat_am$Body_mass_g_RMA)
names(mass_am) <- dat_am$Scientific_name
mass_am <- mass_am[complete.cases(mass_am)]

mass_sq <- log(dat_sq$Body_mass_g_mean)
names(mass_sq) <- dat_sq$Scientific_name
mass_sq <- mass_sq[complete.cases(mass_sq)]

mass_av <- log(dat_av$Body_mass_g_mean)
names(mass_av) <- dat_av$Scientific_name
mass_av <- mass_av[complete.cases(mass_av)]

mass_ma <- log(dat_ma$Body_mass_g_mean)
names(mass_ma) <- dat_ma$Scientific_name
mass_ma <- mass_ma[complete.cases(mass_ma)]

pruned_am <- pruned_sq <- pruned_av <- pruned_ma <- list()

for (i in 1:100) pruned_am[[i]] <-
                    treedata(tr_am[[i]], mass_am, warnings = F)$phy
for (i in 1:100) pruned_sq[[i]] <-
                    treedata(tr_sq[[i]], mass_sq, warnings = F)$phy
for (i in 1:100) pruned_av[[i]] <-
                    treedata(tr_av[[i]], mass_av, warnings = F)$phy
for (i in 1:100) pruned_ma[[i]] <-
                    treedata(tr_ma[[i]], mass_ma, warnings = F)$phy

mass_ord_am <- lapply(pruned_am, function(x) mass_am[x$tip.label])
mass_ord_sq <- lapply(pruned_sq, function(x) mass_sq[x$tip.label])
mass_ord_av <- lapply(pruned_av, function(x) mass_av[x$tip.label])
mass_ord_ma <- lapply(pruned_ma, function(x) mass_ma[x$tip.label])

# function modified from Murrell (2018; DOI: 10.1111/2041-210X.13006)
rank_env_dtt <- function(x, test = "two.sided") {
  
  spp_num <- length(x$times)		
  sims <- x$sim
  sims <- as.matrix(sims)
  
  s1 <- sims[-c(1), ]
  
  r <- x$times[-c(1)]
  
  r <- as.vector(r)
  
  obs <- as.vector(x$dtt)
  obs <- obs[-c(1)]
  
  c1 <- list(r, obs, s1)
  names(c1) = c("r", "obs", "sim_m") 
  c2 <- create_curve_set(c1)
  
  res <- rank_envelope(c2, alternative=test)
  
  final <- list()
  final$times <- res$r
  final$obs <- res$obs
  final$sim <- res$central
  final$rank_lo <- res$lo
  final$rank_hi <- res$hi
  
  return(final)
}

dtt_am <- dtt_sq <- dtt_av <- dtt_ma <- list()

for (i in 1:100) {
  print(i)
  x <- dtt(phy = pruned_am[[i]], data = mass_ord_am[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_am[[i]] <- rank_env_dtt(x)
}
#save(dtt_am, file="data/amphibia/dtt_am_17jan24.RData")

for (i in 1:100) {
  print(i)
  x <- dtt(phy = pruned_sq[[i]], data = mass_ord_sq[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_sq[[i]] <- rank_env_dtt(x)
}
#save(dtt_sq, file="data/reptilia/dtt_sq_17jan24.RData")

for (i in 1:100) {
  print(i)
  x <- dtt(phy = pruned_av[[i]], data = mass_ord_av[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_av[[i]] <- rank_env_dtt(x)
}
#save(dtt_av, file="data/aves/dtt_av_17jan24.RData")

for (i in 1:100) {
  print(i)
  x <- dtt(phy = pruned_ma[[i]], data = mass_ord_ma[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_ma[[i]] <- rank_env_dtt(x)
}
#save(dtt_ma, file="data/mammalia/dtt_ma_17jan24.RData")

#load(file="data/amphibia/dtt_am_17jan24.RData")
#load(file="data/reptilia/dtt_sq_17jan24.RData")
#load(file="data/aves/dtt_av_17jan24.RData")
#load(file="data/mammalia/dtt_ma_17jan24.RData")

cols <- c(rgb(139/255, 71/255, 137/255, 0.4), rgb(97/255, 82/255, 190/255, 0.4),
          rgb(161/255, 204/255, 89/255, 0.5), rgb(255/255, 207/255, 83/255, 0.5))
cols2 <- c(rgb(139/255, 71/255, 137/255), rgb(97/255, 82/255, 190/255),
           rgb(161/255, 204/255, 89/255), rgb(255/255, 207/255, 83/255))
col_sim_pol <- rgb(0.8, 0.8, 0.8, 0.2)
col_sim_cen <- rgb(0.6, 0.6, 0.6, 0.8)

pdf("EvolBiol_20oct23/rev1/Figure3_test.pdf", width = 9.5)

layout(matrix(c(1, 2, 3, 4, 5, 5, 5, 5), ncol = 2, byrow = T),
       heights = c(7, 7, 0.5))

par(mar = c(3, 4, 3, 2))

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "", ylab = "disparity", main = "",
     ylim = c(0, 7))
for (i in 1:100) polygon(c(dtt_am[[i]]$times, rev(dtt_am[[i]]$times)), 
                       c(dtt_am[[i]]$rank_hi, rev(dtt_am[[i]]$rank_lo)), 
                       col = col_sim_pol, border = NA)
for (i in 1:100) lines(dtt_am[[i]]$sim ~ dtt_am[[i]]$times, col = col_sim_cen)
for (i in 1:100) lines(dtt_am[[i]]$obs ~ dtt_am[[i]]$times, col = cols[1])
title("A", adj = 0.95, line = -1.5)

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "", ylab = "", main = "",
     ylim = c(0, 4.1))
for (i in 1:100) polygon(c(dtt_sq[[i]]$times, rev(dtt_sq[[i]]$times)), 
                         c(dtt_sq[[i]]$rank_hi, rev(dtt_sq[[i]]$rank_lo)), 
                         col = col_sim_pol, border = NA)
for (i in 1:100) lines(dtt_sq[[i]]$sim ~ dtt_sq[[i]]$times, col = col_sim_cen)
for (i in 1:100) lines(dtt_sq[[i]]$obs ~ dtt_sq[[i]]$times, col = cols[2])
title("B", adj = 0.95, line = -1.5)

par(mar = c(4, 4, 2, 2))

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", main = "",
     ylab = "disparity", ylim = c(0, 3.5))
for (i in 1:100) polygon(c(dtt_av[[i]]$times, rev(dtt_av[[i]]$times)), 
                         c(dtt_av[[i]]$rank_hi, rev(dtt_av[[i]]$rank_lo)), 
                         col = col_sim_pol, border = NA)
for (i in 1:100) lines(dtt_av[[i]]$sim ~ dtt_av[[i]]$times, col = col_sim_cen)
for (i in 1:100) lines(dtt_av[[i]]$obs ~ dtt_av[[i]]$times, col = cols[3])
title("C", adj = 0.95, line = -1.5)

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", ylab = "",
     main = "", ylim = c(0, 2.5))
for (i in 1:100) polygon(c(dtt_ma[[i]]$times, rev(dtt_ma[[i]]$times)), 
                         c(dtt_ma[[i]]$rank_hi, rev(dtt_ma[[i]]$rank_lo)), 
                         col = col_sim_pol, border = NA)
for (i in 1:100) lines(dtt_ma[[i]]$sim ~ dtt_ma[[i]]$times, col = col_sim_cen)
for (i in 1:100) lines(dtt_ma[[i]]$obs ~ dtt_ma[[i]]$times, col = cols[4])
title("D", adj = 0.95, line = -1.5)

par(mar = c(0, 2, 0, 2))
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1,
     ylim = 0:1)
legend("center", legend  = c("Amphibia", "Squamata", "Aves", "Mammalia"),
       pch = 16, pt.cex = 2, cex = 1.1, bty = "n", horiz  =  T,
       col  =  cols2)

dev.off()
