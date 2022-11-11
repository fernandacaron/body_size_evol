rm(list = ls())

setwd("Documents/Lab/body_size_evol")

library(geiger)

## Carregando os dados
tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")[1:100]
tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")[1:100]
tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")[1:100]
tr_ma <- read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")[1:100]

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

## Fazendo log do tamanho de corpo e pegando so as spp com dados
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

## Tirando as �rvores as spp sem dados
pruned_am <- pruned_sq <- pruned_av <- pruned_ma <- list()

for (i in 1:100) pruned_am[[i]] <-
                    treedata(tr_am[[i]], mass_am, warnings = F)$phy
for (i in 1:100) pruned_sq[[i]] <-
                    treedata(tr_sq[[i]], mass_sq, warnings = F)$phy
for (i in 1:100) pruned_av[[i]] <-
                    treedata(tr_av[[i]], mass_av, warnings = F)$phy
for (i in 1:100) pruned_ma[[i]] <-
                    treedata(tr_ma[[i]], mass_ma, warnings = F)$phy

## Ordenando os dados de tamanho de corpo de acordo com os tips das arvores
mass_ord_am <- lapply(pruned_am, function(x) mass_am[x$tip.label])
mass_ord_sq <- lapply(pruned_sq, function(x) mass_sq[x$tip.label])
mass_ord_av <- lapply(pruned_av, function(x) mass_av[x$tip.label])
mass_ord_ma <- lapply(pruned_ma, function(x) mass_ma[x$tip.label])

dtt_am <- dtt_sq <- dtt_av <- dtt_ma <- list()

dtt_am$dtt <- dtt_am$times <- dtt_am$sim <-
      matrix(nrow = length(mass_am), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = pruned_am[[i]], data = mass_ord_am[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_am$dtt[, i] <- x$dtt
  dtt_am$times[, i] <- x$times
  dtt_am$sim[, i] <- rowMeans(x$sim)
}
#save(dtt_am, file="data/amphibia/dtt_am.RData")

dtt_sq$dtt <- dtt_sq$times <- dtt_sq$sim <-
      matrix(nrow = length(mass_sq), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = pruned_sq[[i]], data = mass_ord_sq[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_sq$dtt[, i] <- x$dtt
  dtt_sq$times[, i] <- x$times
  dtt_sq$sim[, i] <- rowMeans(x$sim)
}
#save(dtt_sq, file="data/reptilia/dtt_sq.RData")

dtt_av$dtt <- dtt_av$times <- dtt_av$sim <-
      matrix(nrow = length(mass_av), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = pruned_av[[i]], data = mass_ord_av[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_av$dtt[, i] <- x$dtt
  dtt_av$times[, i] <- x$times
  dtt_av$sim[, i] <- rowMeans(x$sim)
}
#save(dtt_av, file="data/aves/dtt_av_new.RData")

dtt_ma$dtt <- dtt_ma$times <- dtt_ma$sim <-
      matrix(nrow = length(mass_ma), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = pruned_ma[[i]], data = mass_ord_ma[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  dtt_ma$dtt[, i] <- x$dtt
  dtt_ma$times[, i] <- x$times
  dtt_ma$sim[, i] <- rowMeans(x$sim)
}
#save(dtt_ma, file="data/mammalia/dtt_ma.RData")

#load(file="data/amphibia/dtt_am.RData")
#load(file="data/reptilia/dtt_sq.RData")
#load(file="data/aves/dtt_av_new.RData")
#load(file="data/mammalia/dtt_ma.RData")

cols <- c(rgb(139/255, 71/255, 137/255, 0.4), rgb(97/255, 82/255, 190/255, 0.4),
          rgb(161/255, 204/255, 89/255, 0.5), rgb(255/255, 207/255, 83/255, 0.5))
cols2 <- c(rgb(139/255, 71/255, 137/255), rgb(97/255, 82/255, 190/255),
           rgb(161/255, 204/255, 89/255), rgb(255/255, 207/255, 83/255))
col_sim <- rgb(0.7, 0.7, 0.7, 0.5)

pdf("figures/Figure3.pdf", width = 9.5)

layout(matrix(c(1, 2, 3, 4, 5, 5, 5, 5), ncol = 2, byrow = T),
       heights = c(7, 7, 0.5))

par(mar = c(3, 4, 3, 2))

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "", ylab = "disparity", main = "",
     ylim = c(0, 2))
for (i in 1:100) lines(dtt_am$sim[, i]~dtt_am$times[, i], col = col_sim)
for (i in 1:100) lines(dtt_am$dtt[, i]~dtt_am$times[, i], col = cols[1])
title("A", adj = 0.95, line = -1.5)

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "", ylab = "", main = "",
     ylim = c(0, 2))
for (i in 1:100) lines(dtt_sq$sim[, i]~dtt_sq$times[, i], col = col_sim)
for (i in 1:100) lines(dtt_sq$dtt[, i]~dtt_sq$times[, i], col = cols[2])
title("B", adj = 0.95, line = -1.5)

par(mar = c(4, 4, 2, 2))

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", main = "",
     ylab = "disparity", ylim = c(0, 1.5))
for (i in 1:100) lines(dtt_av$sim[, i]~dtt_av$times[, i], col = col_sim)
for (i in 1:100) lines(dtt_av$dtt[, i]~dtt_av$times[, i], col = cols[3])
title("C", adj = 0.95, line = -1.5)

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", ylab = "",
     main = "", ylim = c(0, 1.5))
for (i in 1:100) lines(dtt_ma$sim[, i]~dtt_ma$times[, i], col = col_sim)
for (i in 1:100) lines(dtt_ma$dtt[, i]~dtt_ma$times[, i], col = cols[4])
title("D", adj = 0.95, line = -1.5)

par(mar = c(0, 2, 0, 2))
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1,
     ylim = 0:1)
legend("center", legend  = c("Amphibia", "Squamata", "Aves", "Mammalia"),
       pch = 16, pt.cex = 2, cex = 1.1, bty = "n", horiz  =  T,
       col  =  cols2)

dev.off()

## SVL

## Fazendo log do SVL e pegando s� as esp�cies com dados
svl_am <- log(dat_am$SVL_mm_1)
names(svl_am) <- dat_am$Scientific_name
svl_am <- svl_am[complete.cases(svl_am)]

svl_sq <- log(dat_sq$SVL_mm_mean)
names(svl_sq) <- dat_sq$Scientific_name
svl_sq <- svl_sq[complete.cases(svl_sq)]

## Tirando as �rvores as spp sem dados
svl_pruned_am <- svl_pruned_sq <- list()

for (i in 1:100) svl_pruned_am[[i]] <-
                    treedata(tr_am[[i]], svl_am, warnings = F)$phy
for (i in 1:100) svl_pruned_sq[[i]] <-
                    treedata(tr_sq[[i]], svl_sq, warnings = F)$phy

## Ordenando os dados de tamanho de corpo de acordo com os tips das �rvores
svl_ord_am <- lapply(svl_pruned_am, function(x) svl_am[x$tip.label])
svl_ord_sq <- lapply(svl_pruned_sq, function(x) svl_sq[x$tip.label])

svl_dtt_am <- svl_dtt_sq <- list()

svl_dtt_am$dtt <- svl_dtt_am$times <- svl_dtt_am$sim <-
      matrix(nrow = length(svl_am), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = svl_pruned_am[[i]], data = svl_ord_am[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  svl_dtt_am$dtt[, i] <- x$dtt
  svl_dtt_am$times[, i] <- x$times
  svl_dtt_am$sim[, i] <- rowMeans(x$sim)
}
#save(svl_dtt_am, file = "data/amphibia/svl_dtt_am.RData")

svl_dtt_sq$dtt <- svl_dtt_sq$times <- svl_dtt_sq$sim <-
      matrix(nrow = length(svl_sq), ncol = 100)
for (i in 1:100) {
  x <- dtt(phy = svl_pruned_sq[[i]], data = svl_ord_sq[[i]], index = "avg.sq",
           nsim = 100, plot = F)
  svl_dtt_sq$dtt[, i] <- x$dtt
  svl_dtt_sq$times[, i] <- x$times
  svl_dtt_sq$sim[, i] <- rowMeans(x$sim)
}
#save(svl_dtt_sq, file = "data/reptilia/svl_dtt_sq.RData")

#load(file = "data/amphibia/svl_dtt_am.RData")
#load(file = "data/reptilia/svl_dtt_sq.RData")

pdf("figures/FigureS6.pdf", width = 9.5)

layout(matrix(c(1, 2, 3, 3), ncol = 2, byrow = T), heights = c(7, 0.5))

par(mar = c(4, 4, 2, 2))

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", main = "",
     ylab = "disparity", ylim = c(0, 1.5))
for (i in 1:100) lines(svl_dtt_am$sim[, i]~svl_dtt_am$times[, i], col = col_sim)
for (i in 1:100) lines(svl_dtt_am$dtt[, i]~svl_dtt_am$times[, i], col = cols[1])
title("A", adj = 0.95, line = -1.5)

plot(0, 1, type = "n", xlim = c(0, 1), xlab = "relative time", ylab = "",
     main = "", ylim = c(0, 1.5))
for (i in 1:100) lines(svl_dtt_sq$sim[, i]~svl_dtt_sq$times[, i], col = col_sim)
for (i in 1:100) lines(svl_dtt_sq$dtt[, i]~svl_dtt_sq$times[, i], col = cols[2])
title("B", adj = 0.95, line = -1.5)

par(mar = c(0, 2, 0, 2))
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1,
     ylim = 0:1)
legend("center", legend  = c("Amphibia", "Squamata"), pch = 16, pt.cex = 2,
       cex = 1.1, bty = "n", horiz = T, col = cols2[1:2])

dev.off()
