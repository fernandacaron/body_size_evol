rm(list = ls())

setwd("Documents/lab/body_size_evol")

#devtools::install_github("mwpennell/arbutus")

library(arbutus)
library(geiger)

## Carregando os dados
tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")
tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_09set21.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

## Fazendo log do tamanho de corpo e pegando só as espécies com dados
mass_am <- log(dat_am$Body_mass_g_1)
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

## Tirando das árvores as spp sem dados 
pruned_am <- pruned_sq <- pruned_av <- pruned_ma <- list()

for (i in 1:1000) pruned_am[[i]] <- treedata(tr_am[[i]], mass_am, warnings = F)$phy
for (i in 1:1000) pruned_sq[[i]] <- treedata(tr_sq[[i]], mass_sq, warnings = F)$phy
for (i in 1:1000) pruned_av[[i]] <- treedata(tr_av[[i]], mass_av, warnings = F)$phy
for (i in 1:1000) pruned_ma[[i]] <- treedata(tr_ma[[i]], mass_ma, warnings = F)$phy

## Ordenando os dados de tamanho de corpo de acordo com os tips das árvores
mass_ord_am <- lapply(pruned_am, function(x) mass_am[x$tip.label])
mass_ord_sq <- lapply(pruned_sq, function(x) mass_sq[x$tip.label])
mass_ord_av <- lapply(pruned_av, function(x) mass_av[x$tip.label])
mass_ord_ma <- lapply(pruned_ma, function(x) mass_ma[x$tip.label])

## Ajustando o modelo de movimento Browniano e vendo adequação do modelo
res_am <- res_sq <- res_av <- res_ma <- matrix(ncol = 13, nrow = 1000)
colnames(res_am) <- colnames(res_sq) <- colnames(res_av) <- colnames(res_ma) <- 
	c("sigsq", 
	  "m.sig.sim", "m.sig.obs", "c.var.sim", "c.var.obs",
      "s.var.sim", "s.var.obs", "s.asr.sim", "s.asr.obs", 
      "s.hgt.sim", "s.hgt.obs", "d.cdf.sim", "d.cdf.obs")

for (i in 1:1000) {
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
#write.csv(res_am, "data/amphibia/res_am.csv")


for (i in 1:1000) {
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
#write.csv(res_av, "data/aves/res_av.csv")


for (i in 1:1000) {
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

#res_am <- read.csv("data/amphibia/res_am.csv")
#res_sq <- read.csv("data/reptilia/res_sq.csv")
#res_av <- read.csv("data/aves/res_av.csv")
#res_ma <- read.csv("data/mammalia/res_ma.csv")

#Figure 4

col_al_sim <- rgb(0.7,0.7,0.7,0.5)
col_sim <- rgb(0.7,0.7,0.7,0.9)
col_al_obs <- c(rgb(139/255,71/255,137/255,0.4), rgb(97/255,82/255,190/255,0.4),
          rgb(161/255,204/255,89/255,0.5), rgb(255/255,207/255,83/255,0.5))
col_obs <- c(rgb(139/255,71/255,137/255), rgb(97/255,82/255,190/255),
          rgb(161/255,204/255,89/255), rgb(255/255,207/255,83/255))

pdf("figures/Figure4.pdf", width = 12, height = 9)

layout(matrix(1:16, ncol = 4, byrow = T))

par(mar = c(2,4,2,2))

#Amphibia
hist(res_am$sigsq, freq = F, xlab = "", col = col_al_obs[1], xlim = c(0, 1),
	 border = col_obs[1], main = "", breaks = 100)

par(mar = c(2,2,2,2))

hist(res_am$c.var.sim, freq = F, xlab = "", ylim = c(0,10), xlim = c(0.65, 3),
     ylab = "", col = col_al_sim, border = col_sim, main = "", breaks = 3)
hist(res_am$c.var.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 100)

hist(res_am$s.asr.sim, freq = F, xlab = "", xlim = c(-0.04, 0.11), ylab = "",
     col = col_al_sim, border = col_sim, main = "", breaks = 20)
hist(res_am$s.asr.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 20)

hist(res_am$s.hgt.sim, freq = F, xlab = "", ylim = c(0,100), 
     xlim = c(-0.2, 0.044), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 20)
hist(res_am$s.hgt.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 20)

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


## SVL

## Fazendo log do tamanho de corpo e pegando só as espécies com dados
svl_am <- log(dat_am$SVL_mm_1)
names(svl_am) <- dat_am$Scientific_name
svl_am <- svl_am[complete.cases(svl_am)]

svl_sq <- log(dat_sq$SVL_mm_mean)
names(svl_sq) <- dat_sq$Scientific_name
svl_sq <- svl_sq[complete.cases(svl_sq)]

## Tirando das árvores as spp sem dados 
svl_pruned_am <- svl_pruned_sq <- list()

for (i in 1:1000) svl_pruned_am[[i]] <- 
  treedata(tr_am[[i]], svl_am, warnings = F)$phy
for (i in 1:1000) svl_pruned_sq[[i]] <- 
  treedata(tr_sq[[i]], svl_sq, warnings = F)$phy

## Ordenando os dados de tamanho de corpo de acordo com os tips das árvores
svl_ord_am <- lapply(svl_pruned_am, function(x) svl_am[x$tip.label])
svl_ord_sq <- lapply(svl_pruned_sq, function(x) svl_sq[x$tip.label])

## Ajustando o modelo de movimento Browniano e vendo adequação do modelo
svl_res_am <- svl_res_sq <- matrix(ncol = 13, nrow = 1000)
colnames(svl_res_am) <- colnames(svl_res_sq) <- 
	c("sigsq", 
	  "m.sig.sim", "m.sig.obs", "c.var.sim", "c.var.obs", 
	  "s.var.sim", "s.var.obs", "s.asr.sim", "s.asr.obs", 
	  "s.hgt.sim", "s.hgt.obs", "d.cdf.sim", "d.cdf.obs")

for (i in 1:1000) {
  svl_fit_bm_am <- fitContinuous(svl_pruned_am[[i]], svl_ord_am[[i]], model = "BM")
  svl.ad.fit_bm_am <- arbutus(svl_fit_bm_am, nsim = 1)
  
  svl_res_am[i, 1] <- svl_fit_bm_am$opt$sigsq
  svl_res_am[i, 2] <- svl.ad.fit_bm_am$sim$m.sig
  svl_res_am[i, 3] <- svl.ad.fit_bm_am$obs$m.sig
  svl_res_am[i, 4] <- svl.ad.fit_bm_am$sim$c.var
  svl_res_am[i, 5] <- svl.ad.fit_bm_am$obs$c.var
  svl_res_am[i, 6] <- svl.ad.fit_bm_am$sim$s.var
  svl_res_am[i, 7] <- svl.ad.fit_bm_am$obs$s.var
  svl_res_am[i, 8] <- svl.ad.fit_bm_am$sim$s.asr
  svl_res_am[i, 9] <- svl.ad.fit_bm_am$obs$s.asr
  svl_res_am[i, 10] <- svl.ad.fit_bm_am$sim$s.hgt
  svl_res_am[i, 11] <- svl.ad.fit_bm_am$obs$s.hgt
  svl_res_am[i, 12] <- svl.ad.fit_bm_am$sim$d.cdf
  svl_res_am[i, 13] <- svl.ad.fit_bm_am$obs$d.cdf
}
#write.csv(svl_res_am, "data/amphibia/svl_res_am.csv")


for (i in 1:1000) {
  svl_fit_bm_sq <- fitContinuous(svl_pruned_sq[[i]], svl_ord_sq[[i]], model = "BM")
  svl.ad.fit_bm_sq <- arbutus(svl_fit_bm_sq, nsim = 1)
  
  svl_res_sq[i, 1] <- svl_fit_bm_sq$opt$sigsq
  svl_res_sq[i, 2] <- svl.ad.fit_bm_sq$sim$m.sig
  svl_res_sq[i, 3] <- svl.ad.fit_bm_sq$obs$m.sig
  svl_res_sq[i, 4] <- svl.ad.fit_bm_sq$sim$c.var
  svl_res_sq[i, 5] <- svl.ad.fit_bm_sq$obs$c.var
  svl_res_sq[i, 6] <- svl.ad.fit_bm_sq$sim$s.var
  svl_res_sq[i, 7] <- svl.ad.fit_bm_sq$obs$s.var
  svl_res_sq[i, 8] <- svl.ad.fit_bm_sq$sim$s.asr
  svl_res_sq[i, 9] <- svl.ad.fit_bm_sq$obs$s.asr
  svl_res_sq[i, 10] <- svl.ad.fit_bm_sq$sim$s.hgt
  svl_res_sq[i, 11] <- svl.ad.fit_bm_sq$obs$s.hgt
  svl_res_sq[i, 12] <- svl.ad.fit_bm_sq$sim$d.cdf
  svl_res_sq[i, 13] <- svl.ad.fit_bm_sq$obs$d.cdf
}
#write.csv(svl_res_sq, "data/reptilia/svl_res_sq.csv")

#svl_res_am <- read.csv("data/amphibia/svl_res_am.csv")
#svl_res_sq <- read.csv("data/reptilia/svl_res_sq.csv")


#########

## SVL

pdf("figures/FigureS7.pdf", width = 12, height = 6)

layout(matrix(1:8, ncol = 4, byrow = T))

par(mar = c(2,4,2,2))

#Amphibia
hist(svl_res_am$sigsq, freq = F, xlab = "", col = col_al_obs[1], xlim = c(0, 1),
	 border = col_obs[1], main = "", breaks = 250)
	 
par(mar = c(2,2,2,2))

hist(svl_res_am$c.var.sim, freq = F, xlab = "", ylim = c(0, 20), 
     xlim = c(0.65, 4), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 1)
hist(svl_res_am$c.var.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1], breaks = 150)
	 
hist(svl_res_am$s.asr.sim, freq = F, xlab = "", ylim = c(0, 100), 
     xlim = c(-0.02, 0.2), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 10)
hist(svl_res_am$s.asr.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1],breaks = 47)

hist(svl_res_am$s.hgt.sim, freq = F, xlab = "", ylim = c(0, 30), 
     xlim = c(-0.68, 0.07), ylab = "", col = col_al_sim, border = col_sim, 
     main = "", breaks = 5)
hist(svl_res_am$s.hgt.obs, freq = F, add = T, col = col_al_obs[1], 
     border = col_obs[1],breaks = 26)

#Squamata
par(mar = c(4,4,2,2))

hist(svl_res_sq$sigsq, freq = F, xlab = expression(sigma^2), xlim = c(0, 1),
     col = col_al_obs[2], border = col_obs[2], main = "", breaks = 3)

par(mar = c(4,2,2,2))

hist(svl_res_sq$c.var.sim, freq = F, xlab = expression("C"["var"]), 
     ylim = c(0, 20), xlim = c(0.65, 4), ylab = "", col = col_al_sim, 
     border = col_sim, main = "", breaks = 1)
hist(svl_res_sq$c.var.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 20)
	 
hist(svl_res_sq$s.asr.sim, freq = F, xlab = expression("S"["asr"]), 
     ylim = c(0, 100), xlim = c(-0.02, 0.2), ylab = "", col = col_al_sim, 
     border = col_sim, main = "", breaks = 5)
hist(svl_res_sq$s.asr.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 20)

hist(svl_res_sq$s.hgt.sim, freq = F, xlab = expression("S"["hgt"]), 
     ylim = c(0, 35), xlim = c(-0.68, 0.07), ylab = "", col = col_al_sim, 
     border = col_sim, main = "", breaks = 5)
hist(svl_res_sq$s.hgt.obs, freq = F, add = T, col = col_al_obs[2], 
     border = col_obs[2], breaks = 26)
	 
dev.off()

