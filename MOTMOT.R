rm(list = ls())

setwd("Documents/lab/body_size_evol")

library(motmot.2.0)
library(phytools)
library(geiger)

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_09set21.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")
num <- sample(1:1000, 10, FALSE)
tr_am <- tr_am[c(num)]
for(i in 1:10) write.tree(tr_am[[i]], "data/amphibia/amphibia_tr10.nex",
                          append = T)

tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")
num <- sample(1:1000, 10, FALSE)
tr_sq <- tr_sq[c(num)]
for(i in 1:10) write.tree(tr_sq[[i]], "data/reptilia/squamata_tr10.nex",
                          append = T)

tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
num <- sample(1:1000, 10, FALSE)
tr_av <- tr_av[c(num)]
for(i in 1:10) write.tree(tr_av[[i]], "data/aves/aves_tr10.nex", append = T)

tr_ma <- read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")
num <- sample(1:1000, 10, FALSE)
tr_ma <- tr_ma[c(num)]
for(i in 1:10) write.tree(tr_ma[[i]], "data/mammalia/mammalia_tr10.nex",
                          append = T)

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

phy_am <- treedata(tr_am[[1]], mass_am, warnings = F)$phy
phy_sq <- treedata(tr_sq[[1]], mass_sq, warnings = F)$phy
phy_av <- treedata(tr_av[[1]], mass_av, warnings = F)$phy
phy_ma <- treedata(tr_ma[[1]], mass_ma, warnings = F)$phy

mass_am <- as.matrix(mass_am[phy_am$tip.label])
mass_sq <- as.matrix(mass_sq[phy_sq$tip.label])
mass_av <- as.matrix(mass_av[phy_av$tip.label])
mass_ma <- as.matrix(mass_ma[phy_ma$tip.label])

tm2_am <- transformPhylo.ML(y = mass_am, phy = phy_am, model = "tm2", 
                            minCladeSize = 10)
tm2_sq <- transformPhylo.ML(y = mass_sq, phy = phy_sq, model = "tm2", 
                            minCladeSize = 10)
tm2_av <- transformPhylo.ML(y = mass_av, phy = phy_av, model = "tm2", 
                            minCladeSize = 10)
tm2_ma <- transformPhylo.ML(y = mass_ma, phy = phy_ma, model = "tm2", 
                            minCladeSize = 10)

summ_am <- traitMedusaSummary(tm2_am, AICc = T)
summ_sq <- traitMedusaSummary(tm2_sq, AICc = T)
summ_av <- traitMedusaSummary(tm2_av, AICc = T)
summ_ma <- traitMedusaSummary(tm2_ma, AICc = T)

pdf("figures/MOTMOT.pdf")

layout(matrix(1:4, ncol = 2), by.row = T)

par(mar = c(1, 1, 1, 1))

plotPhylo.motmot(phy = phy_am, traitMedusaObject = summ_am, reconType = "rates",
                 cex = 0.3, edge.width = 2, show.tip.label = F)
plotPhylo.motmot(phy = phy_sq, traitMedusaObject = summ_sq, reconType = "rates",
                 cex = 0.3, edge.width = 2, show.tip.label = F)
plotPhylo.motmot(phy = phy_av, traitMedusaObject = summ_av, reconType = "rates",
                 cex = 0.3, edge.width = 2, show.tip.label = F)
plotPhylo.motmot(phy = phy_ma, traitMedusaObject = summ_ma, reconType = "rates",
                 cex = 0.3, edge.width = 2, show.tip.label = F)

dev.off()



pdf("figures/MOTMOT_amphibia.pdf")
layout(matrix(1:9, ncol = 3))

par(mar = c(1, 1, 1, 1))

for (i in 2:10) {
    phy <- treedata(tr_am[[i]], mass_am, warnings = F)$phy

    mass_am <- as.matrix(mass_am[phy$tip.label])

    tm2 <- transformPhylo.ML(y = mass_am, phy = phy, model = "tm2", 
                             minCladeSize = 10)

    tm2_summary <- traitMedusaSummary(tm2, AICc = T)

    colour_motmot <- plotPhylo.motmot(phy = phy, 
                                      traitMedusaObject = tm2_summary,
                                      reconType = "rates", cex = 0.3, 
                                      edge.width = 2, show.tip.label = F)
}

dev.off()


