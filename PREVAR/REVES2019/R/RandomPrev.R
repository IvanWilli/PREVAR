# Author: tim
###############################################################################
setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
# ---------------------

lx        <- qx2lx(qx)
# now group / graduate this one
lx_5      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
qx_5      <- lx2qx(lx_5)
prev_5    <- prev_line(lx_5,0,.5)
a5        <- seq(0,110,by=5)
set.seed(1)
LTblock1 <- makeLTblock(lx_5,prev_5,radix=100,interval=5)
LTblock2 <- makeLTblock(lx_5,prev_5,radix=100,interval=5)
LTblock3 <- makeLTblock(lx_5,prev_5,radix=100,interval=5)



pdf("REVES2019/Figs/Random1.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock1,lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()

pdf("REVES2019/Figs/Random2.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock2,lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()

pdf("REVES2019/Figs/Random3.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock3,lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()

# ---------------------------------------
# asymptotic var on .05 grid
lx        <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=.05,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
a5         <- seq(0,110,by=.05)

do.this <- FALSE
if (do.this){
	N <- 1000
	vari <- rep(NA,N)
	for (i in 1:N){
		LTbi <- makeLTblock(lx=lx_05,prev=prev_05,radix=1e5,interval=.05)
		vari[i] <- varblock(LTbi)
		rm(LTbi);gc()
	}
	saveRDS(vari,file="Data/vari1000_1e5_05.rds")
}

#vari <- readRDS("Data/vari1000_1e5_05.rds")
#dput(round(mean(vari),2))
23.35




