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
tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")

tm <- function(taxon, class, y, trait, phy, model, minCladeSize) {
  
  if (class == "Reptilia") {
    subdat <- subset(y, Taxon == taxon, select = trait)[, 1]
    names(subdat) <- subset(y, Taxon == taxon, 
                               select = "Scientific_name")[, 1]
    subdat <- log(subdat)
    subdat <- subdat[complete.cases(subdat)]
  } else {
    subdat <- subset(y, Order == taxon, select = trait)[, 1]
    names(subdat) <- subset(y, Order == taxon, 
                               select = "Scientific_name")[, 1]
    subdat <- log(subdat)
    subdat <- subdat[complete.cases(subdat)]
    }
  
  tr.pruned <- treedata(phy, subdat, warnings = FALSE)$phy

  subdat <- as.matrix(subdat[tr.pruned$tip.label])

  a <- Sys.time()
  tm <- transformPhylo.ML(y = subdat, phy = tr.pruned, model = model,
                           minCladeSize = minCladeSize)
  b <- Sys.time()

  print(b-a) 

  return(tm)

}

plot_tm <- function(summary, col = "blue", cex.label = 0.5, 
                    cex.node = 1.5, edge.width1 = 1, edge.width2 = 2,
                    show.tip.label = TRUE, ...)
{
    nodes <- list()
    for (i in 1:length(as.numeric(summary$Rates[[1]]))) {
        nodes[[i]] <- getDescendants(summary$optimalTree, as.numeric(summary$Rates[[1]][i]))
    }
    nodes <- do.call(c, nodes)
    edge_table <- summary$optimalTree$edge
    rownames(edge_table) <- 1:nrow(edge_table)
    edge_child  <- as.numeric(rownames(edge_table)[edge_table[,2] %in% nodes])

    xxx <- numeric()
    for (i in 1:nrow(edge_table)) xxx[i] <- ifelse(i %in% edge_child, 
                                                   edge.width2, edge.width1)

    plot(summary$optimalTree, show.tip.label = show.tip.label, cex = cex.label,
         edge.width = xxx, ...)
    nodelabels(pch = 16, node = as.numeric(summary$Rates[[1]]), frame = "none",
               col = col, cex = cex.node)
}



tr_am10 <- tm2_anu <- summ_anu <- tm2_cau <- summ_cau <- tm2_gym <- summ_gym <- 
    list()
num <- 1:1000
for (i in 1:1000) {
    err <- try( {
        xxx <- sample(num, 1, FALSE)
        tr_am_test <- tr_am[[xxx]]

        test_tm_anu <- tm(taxon = "Anura", class = "Amphibia", y = dat_am, 
                      trait = "Body_mass_g_1", phy = tr_am_test, 
                      model = "tm2", minCladeSize = 3)
        test_summ_anu <- traitMedusaSummary(test_tm_anu, AICc = T)

        test_tm_cau <- tm(taxon = "Caudata", class = "Amphibia", y = dat_am, 
                      trait = "Body_mass_g_1", phy = tr_am_test, 
                      model = "tm2", minCladeSize = 3)
        test_summ_cau <- traitMedusaSummary(test_tm_cau, AICc = T)

        test_tm_gym <- tm(taxon = "Gymnophiona", class = "Amphibia", 
                          y = dat_am, trait = "Body_mass_g_1", phy = tr_am_test,
                          model = "tm2", minCladeSize = 3)
        test_summ_gym <- traitMedusaSummary(test_tm_gym, AICc = T)
    }
    )

    if (isTRUE(class(err) == "try-error")) { 
        next 
    } else { 
        tr_am10[[length(tr_am10) + 1]] <- tr_am_test

        tm2_anu[[length(tm2_anu) + 1]] <- test_tm_anu
        summ_anu[[length(summ_anu) + 1]] <- test_summ_anu

        tm2_cau[[length(tm2_cau) + 1]] <- test_tm_cau
        summ_cau[[length(summ_cau) + 1]] <- test_summ_cau

        tm2_gym[[length(tm2_gym) + 1]] <- test_tm_gym
        summ_gym[[length(summ_gym) + 1]] <- test_summ_gym
    }

    num <- num[!(num == xxx)]

    if (length(tr_am10) == 10) break
}
#for(i in 1:10) write.tree(tr_am10[[i]], "data/amphibia/amphibia_tr10.tre", 
#                          append = T)
#save(tm2_anu, file = "data/amphibia/tm2_anu.RData")
#save(tm2_cau, file = "data/amphibia/tm2_cau.RData")
#save(tm2_gym, file = "data/amphibia/tm2_gym.RData")

## Gymnophiona: no shifts to plot

pdf("figures/MOTMOT_amphibia.pdf", width = 10)

layout(matrix(1:2, ncol = 2))

par(mar = c(1, 1, 1, 1))

plot_tm(summ_anu[[1]], show.tip.label = T, cex.label = 0.2, edge.width1 = 0.5,
        edge.width2 = 1, cex.node = 1)
title("Anura", cex.main = 1, line = -0.5)
plot_tm(summ_cau[[1]], show.tip.label = T)
title("Caudata", cex.main = 1, line = -0.5)

dev.off()

pdf("figures/MOTMOT_anura.pdf", width = 14, height = 14)

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(0, 0, 0, 0))

plot_tm(summ_anu[[2]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[3]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[4]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[5]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[6]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[7]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[8]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[9]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)
plot_tm(summ_anu[[10]], show.tip.label = T, cex.label = 0.1, edge.width1 = 0.5,
        edge.width2 = 1)

dev.off()

pdf("figures/MOTMOT_caudata.pdf", width = 14, height = 14)

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(1, 1, 1, 1))

plot_tm(summ_cau[[2]], show.tip.label = T)
plot_tm(summ_cau[[3]], show.tip.label = T)
plot_tm(summ_cau[[4]], show.tip.label = T)
plot_tm(summ_cau[[5]], show.tip.label = T)
plot_tm(summ_cau[[6]], show.tip.label = T)
plot_tm(summ_cau[[7]], show.tip.label = T)
plot_tm(summ_cau[[8]], show.tip.label = T)
plot_tm(summ_cau[[9]], show.tip.label = T)
plot_tm(summ_cau[[10]], show.tip.label = T)

dev.off()

########################

## Squamata - Anguimorpha,  Gekkota, Iguania, Lacertoidea, Scincoidea,
## Serpentes

tr_ang10 <- tm2_ang <- summ_ang <- list()
num <- 1:1000
for (i in 1:1000) {
    err <- try( {
        xxx <- sample(num, 1, FALSE)
        tr_sq_test <- tr_sq[[xxx]]

        test_tm <- tm(taxon = "Anguimorpha", class = "Reptilia", y = dat_sq, 
                      trait = "Body_mass_g_mean", phy = tr_sq_test, 
                      model = "tm2", minCladeSize = 3)
        test_summ <- traitMedusaSummary(test_tm, AICc = T)
    }
    )

    if (isTRUE(class(err) == "try-error")) { 
        next 
    } else { 
        tr_ang10[[length(tr_ang10) + 1]] <- tr_sq_test
        tm2_ang[[length(tm2_ang) + 1]] <- test_tm
        summ_ang[[length(summ_ang) + 1]] <- test_summ
    }

    num <- num[!(num == xxx)]

    if (length(tr_ang10) == 10) break
}
#for(i in 1:10) write.tree(tr_ang10[[i]], "data/amphibia/tr_ang10.tre", 
#                          append = T)
#save(tm2_ang, file = "data/reptilia/tm2_ang.RData")

tr_gek10 <- tm2_gek <- summ_gek <- list()
num <- 1:1000
for (i in 1:1000) {
    err <- try( {
        xxx <- sample(num, 1, FALSE)
        tr_sq_test <- tr_sq[[xxx]]

        test_tm <- tm(taxon = "Gekkota", class = "Reptilia", y = dat_sq, 
                      trait = "Body_mass_g_mean", phy = tr_sq_test, 
                      model = "tm2", minCladeSize = 3)
        test_summ <- traitMedusaSummary(test_tm, AICc = T)
    }
    )

    if (isTRUE(class(err) == "try-error")) { 
        next 
    } else { 
        tr_gek10[[length(tr_gek10) + 1]] <- tr_sq_test
        tm2_gek[[length(tm2_gek) + 1]] <- test_tm
        summ_gek[[length(summ_gek) + 1]] <- test_summ
    }

    num <- num[!(num == xxx)]

    if (length(tr_gek10) == 10) break
}
#for(i in 1:10) write.tree(tr_gek10[[i]], "data/amphibia/tr_gek10.tre", 
#                          append = T)
#save(tm2_gek, file = "data/reptilia/tm2_gek.RData")

