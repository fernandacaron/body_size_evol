rm(list = ls())

setwd("Documents/lab/body_size_evol")

library(BAMMtools)
library(coda)

mcmcout <- read.csv("data/aves/bamm/mcmc_out1.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
