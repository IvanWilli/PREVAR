# Author: tim
###############################################################################
setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
# ---------------------
# absurd cases:

lx        <- qx2lx(qx)
# now group / graduate this one
lx_5      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
qx_5      <- lx2qx(lx_5)
prev_5    <- prev_line(lx_5,0,.5)
a5        <- seq(0,110,by=5)

LTblock1 <- LTblock_bottom(lx_5,prev_5,radix=100,interval=5)

pdf("REVES2019/Figs/Bottom1.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock1[100:1,],lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()

###############################################
# Top prev:

LTblock2 <- LTblock_top(lx_5,prev_5,radix=100,interval=5)

pdf("REVES2019/Figs/Top1.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock2[100:1,],lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()

###############################################
# min var

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
prev       <- prev_line(lx,0,.5)
a5         <- seq(0,110,by=5)
a          <- 0:110
LTblock3   <- LTblock_min(lx_05, prev_05,radix=100,interval=5)
#lx=lx_05;prev=prev_05,radix
pdf("REVES2019/Figs/Min1.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_05 / lx_05[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="500 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock3[100:1,],lx=lx_05,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_05,a5,radix=100,lwd=3)
dev.off()



# ----------------------------------------
# vars
lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=.05,method = "monoH.FC"))
prev_05    <- prev_line(lx_05,0,.5)
a5         <- seq(0,110,by=.05)

LTblock1   <- LTblock_bottom(lx_05,prev_05,radix=1e5,interval=.05)
var_max    <- varblock(LTblock1)
rm(LTblock1);gc()

# higher than random
LTblock2   <- LTblock_top(lx_05,prev_05,radix=1e5,interval=.05)
var_top    <- varblock(LTblock2)
rm(LTblock2);gc()

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=.5,method = "monoH.FC"))
prev_05    <- prev_line(lx_05,0,.5)
# lx <-lx_05;prev <- prev_05
a5         <- seq(0,110,by=.5)

LTblock3   <- LTblock_min(lx=lx_05,prev=prev_05,radix=1e3,interval=.5)
hist(rowSums(LTblock3,na.rm=TRUE))
var_min    <- varblock(LTblock3)
rm(LTblock3);gc()
saveRDS(c(top=var_top,max=var_max,min=var_min),file="Data/blockminmax.rds")


# TRUE min should move from omega down, with a target
# almost like a gale-shapely problem?

# radix=200;interval=1;prev <- prev_line(lx,0,.5)
