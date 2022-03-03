rm(list = ls())

setwd("Documents/Lab/body_size_evol")

library(geiger)
library(coda)

## Carregando os dados
tr_am <- read.nexus("data/amphibia/amphibia_VertLife_27JUL20.nex")[1:10]
#tr_sq <- read.nexus("data/reptilia/squamata_VertLife_27JUL20.nex")[1:10]
#tr_av <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")[1:10]
#tr_ma <- 
#  read.nexus("data/mammalia/mammalia_node_dated_VertLife_27JUL20.nex")[1:10]

dat_am <- read.csv("data/amphibia/BodySizeAmphibia_09set21.csv")
#dat_sq <- read.csv("data/reptilia/BodySizeReptilia_15set21.csv")
#dat_av <- read.csv("data/aves/BodySizeAves_09set21.csv")
#dat_ma <- read.csv("data/mammalia/BodySizeMammalia_09set21.csv")

## Selecionando uma árvore aleatória
auteur <- function(tree, data, taxon, col, nchains, ngen, samp, prop.width, 
                   filebase = "result", stats = FALSE) {

  a <- Sys.time()

  data_sep <- subset(data, Order == taxon, select = col)[, 1]
  names(data_sep) <- subset(data, Order == taxon, select = Scientific_name)[, 1]
  data_sep <- data_sep[complete.cases(data_sep)]

  ## Podando árvore 
  pruned <- treedata(tree, data_sep, warnings = F)$phy

  data_order <- data_sep[pruned$tip.label]

  if (missing(prop.width)) {

    prop.width <- calibrate.rjmcmc(pruned, data_order, nstep = 10000, 
                                   type = "jump-rbm")
  
  }

  if (nchains == 2) {
    
    chain1 <- rjmcmc.bm(pruned, data_order, ngen = ngen, samp = samp,
                        prop.width = prop.width, filebase = paste0(filebase, 1),
                        type = "jump-rbm")

    chain2 <- rjmcmc.bm(pruned, data_order, ngen = ngen, samp = samp,
                        prop.width = prop.width, filebase = paste0(filebase, 2),
                        type = "jump-rbm")  

    ps1 <- load.rjmcmc(paste0("jump-relaxedBM.", filebase, 1))

    ps2 <- load.rjmcmc(paste0("jump-relaxedBM.", filebase, 2))
 
    mychains <- mcmc.list(mcmc(ps1$log[,]), mcmc(ps2$log[,])) 

    b <- Sys.time()

    temp <- b - a

    print(temp)

    return(mychains)
 
 } else {

  let <- paste(sample(letters, 9, replace = TRUE), collapse="")

  chain <- rjmcmc.bm(pruned, data_order, ngen = ngen, samp = samp,
                     prop.width = prop.width, filebase = paste0(let, filebase),
                     type = "jump-rbm")

  ps <- load.rjmcmc(paste0("jump-relaxedBM.", paste0(let, filebase)))

  if (stats == TRUE) {

    median_rates <- mean(mcmc(ps$log[, 3]))
    shifts <- mean(mcmc(ps$log[, 4]))
    jumps <- mean(mcmc(ps$log[, 5]))

    b <- Sys.time()

    temp <- b - a 

    print(temp)

    res <- list()

    res$ps <- ps

    res$stats <- data.frame(median_rates = median_rates,
                            shifts = shifts, 
                            jumps = jumps)

    return(res)

  }

    b <- Sys.time()

    temp <- b - a 

    print(temp)

    return(ps)

  }

}

## Amphibia
## Anura(476), Caudata (68), Gymnophiona (8)

### Anura

# Rodando auterur para os anfíbios 5 vezes diferentes com a MESMA árvore para 
# ver se o resultado é consistente independente de desvios

#Time difference of 11.04256 hours
#ngen = 100000000
chains_anu <- auteur(tree = tr_am[[1]], data = dat_am, taxon = "Anura", 
                     col = "Body_mass_g_1", nchains = 2, ngen = 100000000, 
                     samp = 1000, filebase = "anu")

effectiveSize(mcmc(chains_anu[[1]]))
effectiveSize(mcmc(chains_anu[[2]]))


#Potential scale reduction factors:

#        Point est. Upper C.I.
#min           1.00       1.00
#max           1.01       1.01
#median        1.00       1.00
#shifts        1.05       1.22
#jumps         1.00       1.01
#jumpvar       1.01       1.03
#SE            1.01       1.07
#root          1.02       1.05
#lnL           1.01       1.06
#ppos          1.04       1.17

#Multivariate psrf

#1.09

chains_anu10 <- lapply(tr_am[1:10], auteur, data = dat_am, taxon = "Anura", 
                     col = "Body_mass_g_1", nchains = 1, ngen = 50000000, 
                     samp = 1000, filebase = "anus10", stats = TRUE)
res_chain <- as.data.frame(matrix(nrow = 20, ncol = 5))
res_chain[, 1] <- "Amphibia"
res_chain[1:10, 2] <- "Anura"
colnames(res_chain) <- c("Class", "Taxon", "Mean_Median_rates", "Mean_shifts", 
                         "Mean_jumps")
for (i in 1:10) {
  res_chain[i, 3] <- chains_anu10[[i]]$stats[[1]]
  res_chain[i, 4] <- chains_anu10[[i]]$stats[[2]]
  res_chain[i, 5] <- chains_anu10[[i]]$stats[[3]]
}

### Caudata 

chains_cau <- auteur(tree = tr_am[[1]], data = dat_am, taxon = "Caudata", 
                     col = "Body_mass_g_1", nchains = 2, ngen = 100000000, 
                     samp = 1000, filebase = "cau")
effectiveSize(mcmc(chains_cau[[1]]))
effectiveSize(mcmc(chains_cau[[2]]))
plot(mcmc(chains_cau[[1]]))

gelman.diag(chains_cau)

chains_cau10 <- lapply(tr_am[1:10], auteur, data = dat_am, taxon = "Caudata", 
                     col = "Body_mass_g_1", nchains = 1, ngen = 50000000, 
                     samp = 1000, filebase = "caus10", stats = TRUE)

res_chain[11:20, 2] <- "Anura"
for (i in 11:20) {
  res_chain[i, 3] <- chains_cau10[[i]]$stats[[1]]
  res_chain[i, 4] <- chains_cau10[[i]]$stats[[2]]
  res_chain[i, 5] <- chains_cau10[[i]]$stats[[3]]
}


auteur(tr_am, dat_am, "Caudata", "Body_mass_g_1", nloops = 1, ngen = 50000000, 
       samp = 1000, filebase="caudata")
ps_cau1 <- load.rjmcmc("jump-relaxedBM.caudata1")
effectiveSize(ps_cau1$log)
plot(ps_cau1$log, ask = dev.interactive())

ps_anu1 <- load.rjmcmc("jump-relaxedBM.anura1")
ps_anu2 <- load.rjmcmc("jump-relaxedBM.anura2")
ps_anu3 <- load.rjmcmc("jump-relaxedBM.anura3")
ps_anu4 <- load.rjmcmc("jump-relaxedBM.anura4")
ps_anu5 <- load.rjmcmc("jump-relaxedBM.anura5")
effectiveSize(ps_anu1$log)
effectiveSize(ps_anu2$log)
effectiveSize(ps_anu3$log)
effectiveSize(ps_anu4$log)
effectiveSize(ps_anu5$log) # Menor: max (158-349) # Maior jumpvar (8207-9533)

## Squamata

ngen <- 100000000
samp <- 1000

#Time difference of XX
ps_sq1 <- auteur(tr_sq1, mass_sq, ngen = ngen, samp = samp, 
                 filebase="squamata1")
ps_sq2 <- auteur(tr_sq1, mass_sq, ngen = ngen, samp = samp, 
                 filebase="squamata2")
ps_sq3 <- auteur(tr_sq1, mass_sq, ngen = ngen, samp = samp, 
                 filebase="squamata3")
ps_sq4 <- auteur(tr_sq1, mass_sq, ngen = ngen, samp = samp, 
                 filebase="squamata4")
ps_sq5 <- auteur(tr_sq1, mass_sq, ngen = ngen, samp = samp, 
                 filebase="squamata5")

#Time difference of XX
ps_av1 <- auteur(tr_av, mass_av, ngen = ngen, samp = samp, filebase="aves") 
ps_av2 <- auteur(tr_av, mass_av, ngen = ngen, samp = samp, filebase="aves") 
ps_av3 <- auteur(tr_av, mass_av, ngen = ngen, samp = samp, filebase="aves") 
ps_av4 <- auteur(tr_av, mass_av, ngen = ngen, samp = samp, filebase="aves") 
ps_av5 <- auteur(tr_av, mass_av, ngen = ngen, samp = samp, filebase="aves") 

#Time difference of XX
ps_ma1 <- auteur(tr_ma, mass_ma, ngen = ngen, samp = samp, filebase="mammalia") 
ps_ma2 <- auteur(tr_ma, mass_ma, ngen = ngen, samp = samp, filebase="mammalia") 
ps_ma3 <- auteur(tr_ma, mass_ma, ngen = ngen, samp = samp, filebase="mammalia") 
ps_ma4 <- auteur(tr_ma, mass_ma, ngen = ngen, samp = samp, filebase="mammalia") 
ps_ma5 <- auteur(tr_ma, mass_ma, ngen = ngen, samp = samp, filebase="mammalia") 


pdf("figures/auteurAmphibia.pdf")

plot(x = ps_am1, par = "jumps", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am2, par = "jumps", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am3, par = "jumps", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am4, par = "jumps", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am5, par = "jumps", burnin = 0.25, legend = T, show.tip = F)

plot(x = ps_am1, par = "shifts", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am2, par = "shifts", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am3, par = "shifts", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am4, par = "shifts", burnin = 0.25, legend = T, show.tip = F)
plot(x = ps_am5, par = "shifts", burnin = 0.25, legend = T, show.tip = F)

dev.off()

pdf("figures/auteurMCMCAmphibia.pdf")

autocorr.plot(ps_am1$log, ask = dev.interactive())
autocorr.plot(ps_am2$log, ask = dev.interactive())
autocorr.plot(ps_am3$log, ask = dev.interactive())
autocorr.plot(ps_am4$log, ask = dev.interactive())
autocorr.plot(ps_am5$log, ask = dev.interactive())

plot(ps_am1$log, ask = dev.interactive())
plot(ps_am2$log, ask = dev.interactive())
plot(ps_am3$log, ask = dev.interactive())
plot(ps_am4$log, ask = dev.interactive())
plot(ps_am5$log, ask = dev.interactive())

dev.off()

## Squamata
## Anguimorpha (219), Gekkota (1598), Iguania (1753), Lacertoidea  (899),
## Scincoidea (1709), Serpentes (3507)

## Aves
## Apodiformes (435), Charadriiformes (337), Columbiformes (297), 
## Passeriformes (5288), Piciformes (408), Psittaciformes (336)

## Mammalia
## Carnivora (292), Cetartiodactyla (340), Chiroptera (1139), 
## Diprotodontia (146), Eulipotyphla (462), Primates (419), Rodentia (2273)