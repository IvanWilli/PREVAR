
# Author: tim
###############################################################################
setwd("/home/tim/git/PREVAR/PREVAR")

source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")

# try at different granularities:

lx        <- qx2lx(qx)

# now group / graduate this one
lx_5      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
lx_10     <- exp(granularize(log(lx),interval=10,method = "monoH.FC"))
lx.5      <- exp(granularize(log(lx),interval=.5,method = "monoH.FC"))
lx.25     <- exp(granularize(log(lx),interval=.25,method = "monoH.FC"))
lx.1      <- exp(granularize(log(lx),interval=.1,method = "monoH.FC"))
lx.025    <- exp(granularize(log(lx),interval=.025,method = "monoH.FC"))
lx.01     <- exp(granularize(log(lx),interval=.01,method = "monoH.FC"))

# get corresponding qx values
qx_5      <- lx2qx(lx_5)
qx_10     <- lx2qx(lx_10)
qx.5      <- lx2qx(lx.5)
qx.25     <- lx2qx(lx.25)
qx.1      <- lx2qx(lx.1)
qx.025    <- lx2qx(lx.025)
qx.01     <- lx2qx(lx.01)

# and get a prevalence line that goes with each.
prev_10   <- prev_line(lx_10,0,.5)
prev_5    <- prev_line(lx_5,0,.5)
prev      <- prev_line(lx,0,.5)
prev.5    <- prev_line(lx.5,0,.5)
prev.25   <- prev_line(lx.25,0,.5)
prev.1    <- prev_line(lx.1,0,.5)
prev.025  <- prev_line(lx.025,0,.5)
prev.01   <- prev_line(lx.01,0,.5)

# now any comparison can be made with block variance
# versus Bernoulli rewards.










a5      <- seq(0,110,by=5)


LTblock <- makeLTblock(lx_5,prev_5,radix=100,interval=5)
lx_5_100 <- lx_5/(1e5/nrow(LTblock))

plot(a5, lx_5_100, type= 's',
		main = paste0("DLE(0) = ",e0block(LTblock),", sdDLE(0) = ",
				     round(sqrt(varblock(LTblock)),2)))
for (i in 1:nrow(LTblock)){
	plot_traj(LTblock[i,],y=i-1,interval = 5,col=gray(.4))
}
lines(c(0,rep(a5,each=2),110),c(0,rep(prev_5,each=2),.5)*100)
barplot(table(rowSums(LTblock)),space=0)


