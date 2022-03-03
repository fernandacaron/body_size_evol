rm(list = ls())

library(geiger)
library(coda)

## Carregando os dados
tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")[1:10]

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")

subdat <- subset(dat_am, Order == "Anura", select = Body_mass_g_1)[, 1]
names(subdat) <- subset(dat_am, Order == "Anura", select = Scientific_name)[, 1]
subdat <- subdat[complete.cases(subdat)]

pruned <- treedata(tr_am[[1]], subdat, warnings = F)$phy
  
data_order <- subdat[pruned$tip.label]
  
prop.width <- calibrate.rjmcmc(pruned, data_order, nstep = 10000, 
                               type = "jump-rbm")
    
let <- paste(sample(letters, 9, replace = TRUE), collapse="")

chain <- rjmcmc.bm(pruned, data_order, ngen = 100000, samp = 1000,
                   prop.width = prop.width, filebase = paste0(let, "Anu"),
                   type = "jump-rbm")

ps <- load.rjmcmc(paste0("jump-relaxedBM.", paste0(let, "Anu")))
    
effectiveSize(ps$log)
plot(mcmc(ps))
