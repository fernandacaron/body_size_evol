rm(list = ls())

## setwd("Documents/lab/body_size_evol")

library(geiger)
library(coda)

tr <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")[[1]]

dat <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")

subdat <- subset(dat, Order == "Caudata", select = Body_mass_g_1)[, 1]
names(subdat) <- subset(dat, Order == "Caudata", select = Scientific_name)[, 1]
subdat <- subdat[complete.cases(subdat)]

pruned <- treedata(tr, subdat, warnings = F)$phy

subdat <- subdat[pruned$tip.label]

prop.width <- calibrate.rjmcmc(pruned, subdat, nstep = 10000, type = "jump-rbm")

rjmcmc.bm(pruned, subdat, ngen = 5000000, samp = 1000,
          prop.width = prop.width, type = "jump-rbm")

ps <- load.rjmcmc("jump-relaxedBM.result")

effectiveSize(ps$log)
