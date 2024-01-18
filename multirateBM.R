rm(list = ls())

setwd("~/Documents/lab/body_size_evol")

library(phytools)
library(RRphylo)

tree <- pbtree(n = 50)
sig2 <- exp(fastBM(tree, internal = TRUE))
tt <- tree
tt$edge.length <- tt$edge.length*apply(tt$edge, 1, 
                                       function(e,x) mean(x[e]), x = sig2)
x <- fastBM(tt)

multirate <- multirateBM(tree, x, n.iter = 10, parallel = TRUE)
c(RRcova$aces,massdino)->cov.values
c(rownames(RRcova$aces),names(massdino))->names(cov.values)

RRphylo(tree=treedino,y=massdino,cov=cov.values,clus=cc)->RR

RR <- RRphylo(tree = tree, y = x)
cov <- c(RR$aces, x)
names(cov) <- c(rownames(RR$aces), names(x))
RR <- RRphylo(tree = tree, y = x, cov = cov)

rates_bm <- multirate$sig2

rates_rr <- RR$rates[, 1]
names(rates_rr) <- rownames(RR$rates)

rates_bm <- rates_bm[names(sig2)]
rates_rr_ord <- rates_rr[names(sig2)]

sim_sig2 <- apply(tree$edge, 1, function(e, x) mean(x[e]), x = sig2)
mul_sig2 <- apply(tree$edge, 1, function(e, x) mean(x[e]), x = rates_bm)
rrp_sig2 <- apply(tree$edge, 1, function(e, x) mean(x[e]), x = rates_rr_ord)

layout(matrix(1:4, ncol = 2, byrow = T))

par(mar = c(5.1, 5.1, 4.1, 1.1))

plot(sim_sig2, mul_sig2, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     log = "xy", xlim = range(c(sim_sig2, mul_sig2)), 
     ylim = range(c(sim_sig2, mul_sig2)))
lines(range(c(sim_sig2, mul_sig2)), range(c(sim_sig2, mul_sig2)), lwd = 2, 
      col = make.transparent(palette()[4], 0.5))
title(main = "multirateBM (mean)", font.main = 3)
reg_bm1 <- lm(log(mul_sig2) ~ log(sim_sig2))
mtext(paste0("R_sq = ", round(summary(reg_bm1)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)

plot(log(sim_sig2), rrp_sig2, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     xlim = range(c(log(sim_sig2), rrp_sig2)), 
     ylim = range(c(log(sim_sig2), rrp_sig2)))
lines(range(c(log(sim_sig2), rrp_sig2)), range(c(log(sim_sig2), rrp_sig2)), lwd = 2, 
      col = make.transparent(palette()[4], 0.5))
title(main = "RRphylo (mean)", font.main = 3)
reg_rr1 <- lm(rrp_sig2 ~ log(sim_sig2))
mtext(paste0("R_sq = ", round(summary(reg_rr1)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)

sig2_tips <- sig2[names(sig2) %in% tree$tip.label]
rates_bm_tips <- rates_bm[names(rates_bm) %in% tree$tip.label]
rates_bm_tips <- rates_bm_tips[names(sig2_tips)]

plot(sig2_tips, rates_bm_tips, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     log = "xy", xlim = range(c(sig2_tips, rates_bm_tips)), 
     ylim = range(c(sig2_tips, rates_bm_tips)))
lines(range(c(sig2_tips, rates_bm_tips)), range(c(sig2_tips, rates_bm_tips)), 
      lwd = 2, col = make.transparent(palette()[4], 0.5))
title(main = "multirateBM (tips)", font.main = 3)
reg_bm2 <- lm(log(rates_bm_tips) ~ log(sig2_tips))
mtext(paste0("R_sq = ", round(summary(reg_bm2)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)

rates_rr_tips <- rates_rr_ord[names(rates_rr_ord) %in% tree$tip.label]
rates_rr_tips <- rates_rr_tips[names(sig2_tips)]

plot(log(sig2_tips), rates_rr_tips, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     xlim = range(c(log(sig2_tips), rates_rr_tips)), 
     ylim = range(c(log(sig2_tips), rates_rr_tips)))
lines(range(c(log(sig2_tips), rates_rr_tips)), 
      range(c(log(sig2_tips), rates_rr_tips)), lwd = 2, 
      col = make.transparent(palette()[4], 0.5))
title(main = "RRphylo (tips)", font.main = 3)
reg_rr2 <- lm(rates_rr_tips ~ log(sig2_tips))
mtext(paste0("R_sq = ", round(summary(reg_rr2)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)

## with more tips

rm(list = ls())

tree <- pbtree(n = 500)
sig2 <- exp(fastBM(tree, internal = TRUE))
tt <- tree
tt$edge.length <- tt$edge.length*apply(tt$edge, 1, 
                                       function(e,x) mean(x[e]), x = sig2)
x <- fastBM(tt)

trees <- getCladesofSize(tree, 3)

rates_bm <- list()
for (i in 1:length(trees)) {
        print(i)
        dat <- x[names(x) %in% trees[[i]]$tip.label]
        multirate <- multirateBM(trees[[i]], dat, n.iter = 10, parallel = TRUE)

        rates_bm[[i]] <- multirate$sig2
}

rates_bm2 <- do.call(c, rates_bm)
rates_bm2 <- rates_bm2[names(rates_bm2) %in% tree$tip.label]

RR <- RRphylo(tree = tree, y = x)
rates_rr <- RR$rates[, 1]
names(rates_rr) <- rownames(RR$rates)

layout(matrix(1:2, ncol = 2))

par(mar = c(5.1, 5.1, 4.1, 1.1))

sig2_tips <- sig2[names(sig2) %in% tree$tip.label]
rates_bm2 <- rates_bm2[names(sig2_tips)]

plot(sig2_tips, rates_bm2, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     log = "xy", xlim = range(c(sig2_tips, rates_bm2)), 
     ylim = range(c(sig2_tips, rates_bm2)))
lines(range(c(sig2_tips, rates_bm2)), range(c(sig2_tips, rates_bm2)), 
      lwd = 2, col = make.transparent(palette()[4], 0.5))
title(main = "multirateBM (tips)", font.main = 3)
reg_bm2 <- lm(log(rates_bm2) ~ log(sig2_tips))
mtext(paste0("R_sq = ", round(summary(reg_bm2)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)

rates_rr_tips <- rates_rr[names(rates_rr) %in% tree$tip.label]
rates_rr_tips <- rates_rr_tips[names(sig2_tips)]

plot(log(sig2_tips), rates_rr_tips, pch = 21, bg = "grey", cex = 1.5, bty = "n",
     xlab = expression(paste("true rates, ", sigma^2)),
     ylab = expression(paste("estimates rates, ", sigma^2)),
     xlim = range(c(log(sig2_tips), rates_rr_tips)), 
     ylim = range(c(log(sig2_tips), rates_rr_tips)))
lines(range(c(log(sig2_tips), rates_rr_tips)), 
      range(c(log(sig2_tips), rates_rr_tips)), lwd = 2, 
      col = make.transparent(palette()[4], 0.5))
title(main = "RRphylo (tips)", font.main = 3)
reg_rr2 <- lm(rates_rr_tips ~ log(sig2_tips))
mtext(paste0("R_sq = ", round(summary(reg_rr2)$r.squared, 3)), line = -2, 
      side = 3, adj = 0.1, cex = 0.7)







H <- max(nodeHeights(pruned_tr_am))
xxx <- treeSlice(tr_slice, slice, trivial = T)

                str(xxx)

                i_start = i
                names_tr <- character()
                for (j in 1:length(xxx)) {
                        if (length(xxx[[j]]$tip.label) <= 50) {
                                tr[[i]] = xxx[[j]]
                                names_tr <- c(names_tr, xxx[[j]]$tip.label)
                                i = i + 1
                        }
                }
                i_end = i

                tr_slice <- drop.tip(tr_slice, names_tr)

        }


        for (i in 1:length(tr_slice))
        xxx <- treeSlice(tr_slice, slice, trivial = T)

        str(xxx)

        i_start = i
        names_tr <- character()
        for (j in 1:length(xxx)) {
                if (length(xxx[[j]]$tip.label) <= 50) {
                        tr[[i]] = xxx[[j]]
                        names_tr <- c(names_tr, xxx[[j]]$tip.label)
                        i = i + 1
                }
        }
        i_end = i

        tr_slice <- drop.tip(tr_slice, names_tr)

        tr_slice

        if (length(tr_slice$tip.label) > 50) {
                continue = TRUE
        } else {
                continue = FALSE
        }

        i
}

xxx2 <- treeSlice(xxx[[1]], 1, trivial=TRUE)
