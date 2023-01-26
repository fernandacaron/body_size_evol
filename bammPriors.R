rm(list = ls())

setwd("Documents/lab/body_size_evol")

library(phytools)
library(BAMMtools)
library(geiger)

#tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")
#tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
#tr_ma <- read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")

#dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
#dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_10set22.csv")
#dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

#mass_am <- dat_am$Body_mass_g_1
#names(mass_am) <- dat_am$Scientific_name

#mass_sq <- dat_sq$Body_mass_g_mean
#names(mass_sq) <- dat_sq$Scientific_name

mass_av <- dat_av$Body_mass_g_mean
names(mass_av) <- dat_av$Scientific_name

#mass_ma <- dat_ma$Body_mass_g_mean
#names(mass_ma) <- dat_ma$Scientific_name

#mass_am <- mass_am[complete.cases(mass_am)]
#mass_sq <- mass_sq[complete.cases(mass_sq)]
mass_av <- mass_av[complete.cases(mass_av)]
#mass_ma <- mass_ma[complete.cases(mass_ma)]

for (i in 1:1000) {
  #tr_am[[i]] <- treedata(tr_am[[1]], mass_am)$phy
  #tr_sq[[i]] <- treedata(tr_sq[[1]], mass_sq)$phy
  tr_av[[i]] <- treedata(tr_av[[1]], mass_av)$phy
  #tr_ma[[i]] <- treedata(tr_ma[[1]], mass_ma)$phy
}

#mass_am <- mass_am[tr_am[[1]]$tip.label]
#mass_sq <- mass_sq[tr_sq[[1]]$tip.label]
mass_av <- mass_av[tr_av[[1]]$tip.label]
#mass_ma <- mass_ma[tr_ma[[1]]$tip.label]

for (i in 1:1) {
  priors <- setBAMMpriors(phy = tr_av[[i]], traits = mass_av, outfile = NULL)
  write.tree(tr_av[[i]], paste0("data/aves/bamm/tree", i, ".tre"))
  generateControlFile(file = paste0("data/aves/bamm/control", i, ".txt"),
                      type = "trait",
                      params = list(
                        treefile = paste0("tree", i, ".tre"),
                        traitfile = "trait.txt",
                        runInfoFilename = paste0("run_info", i, ".txt"),
                        eventDataInfile = paste0("event_data_in", i, ".txt"),
                        mcmcOutfile = paste0("mcmc_out", i, ".txt"), 
                        eventDataOutfile = paste0("event_data", i, ".txt"), 
                        chainSwapFileName = paste0("chain_swap", i, ".txt"),
                        numberOfGenerations = "100000000",
                        mcmcWriteFreq = 10000,
                        eventDataWriteFreq = 10000,
                        printFreq = 10000,
                        overwrite = "0",
                        betaInitPrior = as.numeric(priors["betaInitPrior"]),
                        betaShiftPrior = as.numeric(priors["betaShiftPrior"]),
                        useObservedMinMaxAsTraitPriors = as.numeric(priors["useObservedMinMaxAsTraitPriors"]),
                        expectedNumberOfShifts = "50"))
}

