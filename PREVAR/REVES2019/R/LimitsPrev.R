# Author: tim
###############################################################################
setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
# ---------------------
# absurd cases:
LTblock_bottom <- function(lx, prev, radix = 1e5, interval = 1){
	lx    <- round(lx / lx[1] * radix)
	
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		Pmat[(radix-Nsick[i]):radix,i] <- interval
	}
	Pmat
}

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
LTblock_top <- function(lx, prev, radix = 1e5, interval = 1){
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		if (Nsick[i] > 0){
			Pmat[(ind+1):(ind+Nsick[i]),i] <- interval
		}
	}
	Pmat
}

LTblock2 <- LTblock_top(lx_5,prev_5,radix=100,interval=5)

pdf("REVES2019/Figs/Top1.pdf")
par(xaxs="i",yaxs="i")
plot(a5, round(lx_5 / lx_5[1] * 100), type = 's', las = 1,
		xlab="Age",ylab="100 lives",cex.lab=1.5,cex.axis=1.5,lwd=2,ylim=c(0,101))
plot_LTblock(LTblock2[100:1,],lx=lx_5,a=a5,interval=5,add=TRUE,col="black")
lx_block_line(lx_5,a5,radix=100,lwd=3)
dev.off()


# ----------------------------------------
# vars
lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=.05,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
a5         <- seq(0,110,by=.05)

LTblock1   <- LTblock_bottom(lx_05,prev_05,radix=1e5,interval=.05)
var_max    <- varblock(LTblock1)
rm(LTblock1);gc()

# WHY higher than random?
LTblock2   <- LTblock_top(lx_05,prev_05,radix=1e5,interval=.05)
var_min    <- varblock(LTblock2)
rm(LTblock2);gc()

LTbi <- makeLTblock(lx=lx_05,prev=prev_05,radix=1e5,interval=.05)
varblock(LTbi)
rm(LTbi);gc()
plot(sort(rowSums(LTbi)))
lines(sort(rowSums(LTblock2,na.rm=TRUE)))


# TRUE min should move from omega down, with a target
# almost like a gale-shapely problem?

# radix=200;interval=1;prev <- prev_line(lx,0,.5)
LTblock_min <- function(lx, prev, radix = 100, interval = 1){
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		if (Nsick[i] > 0){
			Pmat[(ind+1):(ind+Nsick[i]),i] <- interval
		}
	}
	
	# now for the shuffle, code not optimal
	var_i           <- varblock(Pmat)
	for (i in 1:(sum(Nsick)/2)){
		var_old         <- var_i
		Di              <- rowSums(Pmat, na.rm = TRUE)
		resids          <- Di - mean(Di)
		posi            <- resids > 0
		negi            <- resids < 0
		prob_from       <- resids
		prob_to         <- resids
		prob_from[negi] <- 0
		
		# only go down
		prob_to[posi]   <- 0
		prob_to[1:from_i]  <- 0
		prob_to         <- abs(prob_to)
		
		# picks out cells on bottom, when dropping from top that are 
		# also in the leftmost position
		pd              <- diff(rbind(Pmat,interval))
		# -interval happens when dropping top down
		from_cells      <- pd == -interval 
		from_i          <- which.max(prob_from)
		from_col        <- which(diff(c(0,from_cells[from_i,])) == interval)
		Pmat[from_i,]
		prob_to[Pmat[,from_col] == interval & !is.na(Pmat[,from_col])] <- 0
		Pmat[from_i,from_col] <- 0
		Pmat[which.max(prob_to),from_col] <- interval
		var_i           <- varblock(Pmat)
		if (!(var_i <= var_old)){
			break
		}
	}
	
}