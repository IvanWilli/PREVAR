
# Author: tim
###############################################################################

setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
source("R/BlockFunctions.R")
#----------------------

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx), interval = .05, method = "monoH.FC"))  
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05, 0, .5)

var_fixed <- VarFixed(qx_05,prev_05,interval=.05,calcs = "faster",closeout=TRUE)
var_bern  <- VarBernoulli(qx_05,prev_05,interval=.05,calcs = "faster",closeout=TRUE)

saveRDS(c(bern=var_bern,fixed=var_fixed),"Data/varbernfixed.rds")

###########################################################
# need code to draw fixed rewards prevalence, same 5-year 100 radix setup.

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
a5         <- seq(0,110,by=5)

cols       <- gray(seq(1,0,length=length(a5)))
lx_05      <- round(lx_05 / lx_05[1] * 100)
#lx=lx_05;prev=prev_05,radix
pdf("REVES2019/Figs/FixedDecimal.pdf")
par(xaxs="i",yaxs="i")
plot(a5, lx_05, type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
for (i in 1:(length(a5)-1)){
	rect(a5[i],0,a5[i+1],lx_05[i],col=cols[i],border=NA)
}
lx_block_line(lx_05,a5,radix=100,lwd=3)
dev.off()
