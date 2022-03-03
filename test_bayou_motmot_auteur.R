## Testando bayou

rm(list =l s())

setwd("Documents/Lab/body_size_evol")

library(bayou)

tr <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")[[1]]

dat <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")

mass <- log(dat$Body_mass_g_mean)
names(mass) <- dat$Scientific_name
mass <- mass[complete.cases(mass)]

pruned.tr <- treedata(tr, mass, warnings = F)$phy

mass <- mass[pruned.tr$tip.label]

## Prior: is going to take our parameters and output the prior probability of 
## our parameter values. It represents the initial degree of belief in what 
## values the parameters will take
## One trick to make sure the prior functions are reasonable is to simulate a 
## bunch of values and take the quantiles of the distribution

priorOU <- make.prior(pruned.tr) #default values

## Initiate the MCMC chain with some starting values
## It's good to run multiple chains from multiple different starting values
## Let's simulate some values from the prior distribution and make sure the 
## priorfunctions works

startpars <- priorSim(priorOU, pruned.tr, plot = TRUE)$pars[[1]]

priorOU(startpars)

## Making an object to manage our MCMC analysis

mcmcOU <- bayou.makeMCMC(pruned.tr, mass, SE = 0.1, prior = priorOU,  
                         outname = "modelOU", plot.freq = NULL) 

## Run the MCMC

mcmcOU$run(10000) 

## The full MCMC results are written to a set of files. Load them back in to R 
## as follows

chainOU <- mcmcOU$load()

## Set a "burnin" parameter that tells the package coda to discard the first 
## bit of the chain

chainOU <- set.burnin(chainOU, 0.3)

summary(chainOU)

plot(chainOU, auto.layout = FALSE)

## Note the effective sample sizes in our summary (the NA's for the all theta 
## row are expected, this is because these aren't a single parameter, but a 
## variable number of optima that are coming in and out of existence throughout 
##the chain)

## 3 alternative ways of visualizing our chain

par(mfrow = c(2,2))

par(mar = c(0,0,0,0))

plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3, cex = 0.3)

plotBranchHeatMap(pruned.tr, chainOU, "theta", burnin = 0.3, pal = cm.colors)

phenogram.density(pruned.tr, mass, burnin = 0.3, chainOU, pp.cutoff = 0.3)

#Even though we haven't gotten convergence yet, we're probably picking up the 
#major shifts pretty well.

##################

## Testando MOTMOT

rm(list = ls())

setwd("C:/Users/ferna/Documents/IC/Body_size_evol")

library(motmot)
library(phytools)
library(geiger)

tr <- pbtree(n = 100)
dat <- as.matrix(fastBM(tr))

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
dat_av <- read.csv("data/aves/BodySizeAves_09set21.csv")
dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")
tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexues("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")

tr_1 <- tr_sq[[1]]
mass_sq <- log(dat_sq$Body_mass_g_mean)
names(mass_sq) <- dat_sq$Scientific_name
mass_sq <- mass_sq[complete.cases(mass_sq)]

tr_1 <- treedata(tr_1, mass_sq, warnings = F)$phy

mass_sq <- as.matrix(mass_sq[tr_1$tip.label])

## tm1 algorithm: minCladeSize defines the minimum size of clades on which to 
## infer new rate shifts; nSplits defines the maximum number of rate shifts
a <- Sys.time()
tm1_ml <- transformPhylo.ML(y = mass_sq, phy = tr_1, model = "tm1", 
                            minCladeSize = 10)
b <- Sys.time()
b-a

#minCladeSize=10 Time difference of 3.762855 mins
#minCladeSize=10 Time difference of 7.139096 mins

tm1_summary <- summary.traitMedusa(tm1_ml, cutoff = 2, AICc = T) 

par(mar = c(0, 0, 0, 0))
colour_motmot <- plot.traitMedusa.model(x = tm1_summary, reconType = "rates", 
                                        cex = 0.3)

## tm2 model
a <- Sys.time()
tm2_ml <- transformPhylo.ML(y = mass_am, phy = tr_1, model = "tm2", 
                            minCladeSize = 5)
b <- Sys.time()
b-a

#minCladeSize=10 Time difference of 5.594881 mins
#minCladeSize=5 Time difference of 11.22721 mins

tm2_summary <- traitMedusaSummary(tm2_ml, cutoff = 2, AICc = T)

colour_motmot <- plotPhylo.motmot(phy = tr, traitMedusaObject = tm2_summary, 
                                  reconType = "rates", cex = 0.3, 
                                  edge.width = 2)



##################

## Testando AUTEUR

Examples

library(geiger)

## GENERATE DATA: jump-diffusion
phy <- ladderize(sim.bdtree(n = 200), right = FALSE)
r <- paste(sample(letters,9,replace = TRUE),collapse = "")
defpar <- par(no.readonly = TRUE)


tmp <- ex.jumpsimulator(phy, jumps = 10)
dat <- tmp$dat
hist <- tmp$hist

ex.traitgram(phy, hist, alpha = 0) # plot history of trait change

## RUN ANALYSIS


## coda package is not a dependency of geiger
## but is very useful for evaluating mcmc runs
library(coda)


rjmcmc.bm(phy, dat, prop.width = 1.5, ngen = 20000, samp = 500, filebase = r,
          simple.start = TRUE, type = "jump-bm")
outdir <- paste("jump-BM", r, sep = ".")
ps <- load.rjmcmc(outdir)

dev.new()
plot(x = ps, par = "jumps", burnin = 0.25, legend = FALSE, show.tip = FALSE, type = "fan", edge.width = 2)
mm = match(phy$edge[,2],hist$descendant)
hist = hist[mm,]
edgelabels.auteur(text = NULL, pch = 21, cex = hist$cex, bg = NA, 
                  col = ifelse(hist$cex>0, 1, NA), lty = 2)
title("red (estimated); black (true jump size)", line = -5)
par(defpar)

dev.new()
## from the coda package
coda::autocorr.plot(ps$log, ask = dev.interactive())
plot(ps$log, ask = dev.interactive())

## GENERATE DATA: multi-rate diffusion
scl <- ex.ratesimulator(phy, min = 12, show.tip = FALSE)
dat <- rTraitCont(scl)

## RUN ANALYSIS
rjmcmc.bm(phy, dat, prop.width = 1.5, ngen = 20000, samp = 500, filebase = r, simple.start = TRUE, type = "rbm")
outdir <- paste("relaxedBM", r, sep = ".")
ps <- load.rjmcmc(outdir)
dev.new()
plot(x = ps, par = "shifts", burnin = 0.25, legend = TRUE, show.tip = FALSE, edge.width = 2)

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
unlink(dir(pattern = paste(r)),recursive = TRUE)


