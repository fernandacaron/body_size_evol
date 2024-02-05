rm(list = ls())

setwd("~/Documents/lab/body_size_evol")

#devtools::install_github("mwpennell/arbutus")

library(arbutus)
library(geiger)

tr_am <- read.nexus("~/Documents/lab/data/trees/amphibia_VertLife_27JUL20.nex")
tr_sq <- read.nexus("~/Documents/lab/data/trees/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("~/Documents/lab/data/trees/mammalia_node_dated_VertLife_27JUL20.nex")

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

for (i in 1:1000) pruned_am[[i]] <- treedata(tr_am[[i]], mass_am, warnings = F)$phy
for (i in 1:1000) pruned_sq[[i]] <- treedata(tr_sq[[i]], mass_sq, warnings = F)$phy
for (i in 1:1000) pruned_av[[i]] <- treedata(tr_av[[i]], mass_av, warnings = F)$phy
for (i in 1:1000) pruned_ma[[i]] <- treedata(tr_ma[[i]], mass_ma, warnings = F)$phy

mass_ord_am <- lapply(pruned_am, function(x) mass_am[x$tip.label])
mass_ord_sq <- lapply(pruned_sq, function(x) mass_sq[x$tip.label])
mass_ord_av <- lapply(pruned_av, function(x) mass_av[x$tip.label])
mass_ord_ma <- lapply(pruned_ma, function(x) mass_ma[x$tip.label])

res_am <- res_sq <- res_av <- res_ma <- matrix(ncol = 13, nrow = 1000)
colnames(res_am) <- colnames(res_sq) <- colnames(res_av) <- colnames(res_ma) <- 
	c("sigsq", 
	  "m.sig.sim", "m.sig.obs", "c.var.sim", "c.var.obs",
      "s.var.sim", "s.var.obs", "s.asr.sim", "s.asr.obs", 
      "s.hgt.sim", "s.hgt.obs", "d.cdf.sim", "d.cdf.obs")

for (i in 1:1000) {
     print(i)
  fit_bm_am <- fitContinuous(pruned_am[[i]], mass_ord_am[[i]], model = "BM")
  ad.fit_bm_am <- arbutus(fit_bm_am, nsim = 1)
  
  res_am[i, 1] <- fit_bm_am$opt$sigsq
  res_am[i, 2] <- ad.fit_bm_am$sim$m.sig
  res_am[i, 3] <- ad.fit_bm_am$obs$m.sig
  res_am[i, 4] <- ad.fit_bm_am$sim$c.var
  res_am[i, 5] <- ad.fit_bm_am$obs$c.var
  res_am[i, 6] <- ad.fit_bm_am$sim$s.var
  res_am[i, 7] <- ad.fit_bm_am$obs$s.var
  res_am[i, 8] <- ad.fit_bm_am$sim$s.asr
  res_am[i, 9] <- ad.fit_bm_am$obs$s.asr
  res_am[i, 10] <- ad.fit_bm_am$sim$s.hgt
  res_am[i, 11] <- ad.fit_bm_am$obs$s.hgt
  res_am[i, 12] <- ad.fit_bm_am$sim$d.cdf
  res_am[i, 13] <- ad.fit_bm_am$obs$d.cdf
}
#write.csv(res_am, "data/amphibia/res_am_17jan24.csv")

for (i in 1:1000) {
  print(i)
  fit_bm_sq <- fitContinuous(pruned_sq[[i]], mass_ord_sq[[i]], model = "BM")
  ad.fit_bm_sq <- arbutus(fit_bm_sq, nsim = 1)
  
  res_sq[i, 1] <- fit_bm_sq$opt$sigsq
  res_sq[i, 2] <- ad.fit_bm_sq$sim$m.sig
  res_sq[i, 3] <- ad.fit_bm_sq$obs$m.sig
  res_sq[i, 4] <- ad.fit_bm_sq$sim$c.var
  res_sq[i, 5] <- ad.fit_bm_sq$obs$c.var
  res_sq[i, 6] <- ad.fit_bm_sq$sim$s.var
  res_sq[i, 7] <- ad.fit_bm_sq$obs$s.var
  res_sq[i, 8] <- ad.fit_bm_sq$sim$s.asr
  res_sq[i, 9] <- ad.fit_bm_sq$obs$s.asr
  res_sq[i, 10] <- ad.fit_bm_sq$sim$s.hgt
  res_sq[i, 11] <- ad.fit_bm_sq$obs$s.hgt
  res_sq[i, 12] <- ad.fit_bm_sq$sim$d.cdf
  res_sq[i, 13] <- ad.fit_bm_sq$obs$d.cdf
}
#write.csv(res_sq, "data/reptilia/res_sq.csv")


for (i in 1:1000) {
     print(i)
  fit_bm_av <- fitContinuous(pruned_av[[i]], mass_ord_av[[i]], model = "BM")
  ad.fit_bm_av <- arbutus(fit_bm_av, nsim = 1)
  
  res_av[i, 1] <- fit_bm_av$opt$sigsq
  res_av[i, 2] <- ad.fit_bm_av$sim$m.sig
  res_av[i, 3] <- ad.fit_bm_av$obs$m.sig
  res_av[i, 4] <- ad.fit_bm_av$sim$c.var
  res_av[i, 5] <- ad.fit_bm_av$obs$c.var
  res_av[i, 6] <- ad.fit_bm_av$sim$s.var
  res_av[i, 7] <- ad.fit_bm_av$obs$s.var
  res_av[i, 8] <- ad.fit_bm_av$sim$s.asr
  res_av[i, 9] <- ad.fit_bm_av$obs$s.asr
  res_av[i, 10] <- ad.fit_bm_av$sim$s.hgt
  res_av[i, 11] <- ad.fit_bm_av$obs$s.hgt
  res_av[i, 12] <- ad.fit_bm_av$sim$d.cdf
  res_av[i, 13] <- ad.fit_bm_av$obs$d.cdf
}
#write.csv(res_av, "data/aves/res_av_new.csv")


for (i in 1:1000) {
  print(i)
  fit_bm_ma <- fitContinuous(pruned_ma[[i]], mass_ord_ma[[i]], model = "BM")
  ad.fit_bm_ma <- arbutus(fit_bm_ma, nsim = 1)
  
  res_ma[i, 1] <- fit_bm_ma$opt$sigsq
  res_ma[i, 2] <- ad.fit_bm_ma$sim$m.sig
  res_ma[i, 3] <- ad.fit_bm_ma$obs$m.sig
  res_ma[i, 4] <- ad.fit_bm_ma$sim$c.var
  res_ma[i, 5] <- ad.fit_bm_ma$obs$c.var
  res_ma[i, 6] <- ad.fit_bm_ma$sim$s.var
  res_ma[i, 7] <- ad.fit_bm_ma$obs$s.var
  res_ma[i, 8] <- ad.fit_bm_ma$sim$s.asr
  res_ma[i, 9] <- ad.fit_bm_ma$obs$s.asr
  res_ma[i, 10] <- ad.fit_bm_ma$sim$s.hgt
  res_ma[i, 11] <- ad.fit_bm_ma$obs$s.hgt
  res_ma[i, 12] <- ad.fit_bm_ma$sim$d.cdf
  res_ma[i, 13] <- ad.fit_bm_ma$obs$d.cdf
}
#write.csv(res_ma, "data/mammalia/res_ma.csv")

#res_am <- read.csv("data/amphibia/res_am_17jan24.csv")
#res_sq <- read.csv("data/reptilia/res_sq.csv")
#res_av <- read.csv("data/aves/res_av_new.csv")
#res_ma <- read.csv("data/mammalia/res_ma.csv")

#Figure 4

col_al_sim <- rgb(0.7,0.7,0.7,0.5)
col_sim <- rgb(0.7,0.7,0.7,0.9)
col_al_obs <- c(rgb(139/255,71/255,137/255,0.4), rgb(97/255,82/255,190/255,0.4),
          rgb(161/255,204/255,89/255,0.5), rgb(255/255,207/255,83/255,0.5))
col_obs <- c(rgb(139/255,71/255,137/255), rgb(97/255,82/255,190/255),
          rgb(161/255,204/255,89/255), rgb(255/255,207/255,83/255))

pdf("EvolBiol_20oct23/rev1/Figure4.pdf", width = 12, height = 9)

layout(matrix(1:16, ncol = 4, byrow = T))

par(mar = c(2,4,2,2))

#Amphibia
hist(res_am$sigsq, freq = F, xlab = "", col = col_al_obs[1], xlim = c(0, 1),
	 border = col_obs[1], main = "", breaks = 600)

par(mar = c(2,2,2,2))

hist(res_am$c.var.sim, freq = F, xlab = "", ylim = c(0,15), xlim = c(0.65, 3),
     ylab = "", col = col_al_sim, border = col_sim, main = "", breaks = 1)
hist(res_am$c.var.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 400)

hist(res_am$s.asr.sim, freq = F, xlab = "", xlim = c(-0.13, 0.11), ylab = "",
     col = col_al_sim, border = col_sim, main = "", breaks = 5)
hist(res_am$s.asr.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 30)

hist(res_am$s.hgt.sim, freq = F, xlab = "", ylim = c(0,50), 
     xlim = c(-0.3, 0.044), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 20)
hist(res_am$s.hgt.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 60)

#Squamata
par(mar = c(2,4,2,2))

hist(res_sq$sigsq, freq = F, xlab = "", col = col_al_obs[2], xlim = c(0, 1),
     border = col_obs[2], main = "", breaks = 20)

par(mar = c(2,2,2,2))
	 
hist(res_sq$c.var.sim, freq = F, xlab = "", xlim = c(0.65, 3),
     ylab = "", col = col_al_sim, border = col_sim, main = "", breaks = 1)
hist(res_sq$c.var.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 20)

hist(res_sq$s.asr.sim, freq = F, xlab = "", xlim = c(-0.04, 0.11), 
     ylab = "", col = col_al_sim, border = col_sim, main = "", breaks = 2)
hist(res_sq$s.asr.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 10)

hist(res_sq$s.hgt.sim, freq = F, xlab = "", ylim = c(0, 150), 
     xlim = c(-0.2, 0.044), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 3)
hist(res_sq$s.hgt.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 15)

#Aves
par(mar = c(2,4,2,2))

hist(res_av$sigsq, freq = F, xlab = "", col = col_al_obs[3], xlim = c(0, 1),
     border = col_obs[3], main = "", breaks = 350)

par(mar = c(2,2,2,2))

hist(res_av$c.var.sim, freq = F, xlab = "", xlim = c(0.65, 3), ylab = "", 
     col = col_al_sim, border = col_sim, main = "", breaks = 1)
hist(res_av$c.var.obs, freq = F, add = T, col = col_al_obs[3], 
     border = col_obs[3], breaks = 400)
	 
hist(res_av$s.asr.sim, freq = F, xlab = "", xlim = c(-0.04, 0.11), ylab = "",
     col = col_al_sim, border = col_sim, main = "", breaks = 5)
hist(res_av$s.asr.obs, freq = F, add = T, col = col_al_obs[3], 
     border = col_obs[3], breaks = 15)

hist(res_av$s.hgt.sim, freq = F, xlab = "", ylim = c(0, 120), 
     xlim = c(-0.2, 0.044), col = col_al_sim, border = col_sim, main = "",
     breaks = 10)
hist(res_av$s.hgt.obs, freq = F, add = T, col = col_al_obs[3], 
     border = col_obs[3], breaks = 30)
	 
#Mammalia
par(mar = c(4,4,2,2))

hist(res_ma$sigsq, freq = F, xlab = expression(sigma^2), col = col_al_obs[4], 
     border = col_obs[4], main = "", breaks = 10, xlim = c(0, 1))

par(mar = c(4,2,2,2))

hist(res_ma$c.var.sim, freq = F, xlab = expression("C"["var"]), 
     xlim = c(0.65, 3), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 1)
hist(res_ma$c.var.obs, freq = F, add = T, col = col_al_obs[4], 
     border = col_obs[4], breaks = 20)
	 
hist(res_ma$s.asr.sim, freq = F, xlab = expression("S"["asr"]), 
     xlim = c(-0.04, 0.11), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 5)
hist(res_ma$s.asr.obs, freq = F, add = T, col = col_al_obs[4], 
     border = col_obs[4], breaks = 5)

hist(res_ma$s.hgt.sim, freq = F, xlab = expression("S"["hgt"]), 
     ylim = c(0, 120), xlim = c(-0.2, 0.044), ylab = "", col = col_al_sim,
     border = col_sim, main = "", breaks = 10)
hist(res_ma$s.hgt.obs, freq = F, add = T, col = col_al_obs[4], 
     border = col_obs[4], breaks = 15)

dev.off()

