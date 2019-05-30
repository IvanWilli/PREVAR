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

###############################################
# min var
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
	#var_i           <- varblock(Pmat)
	NN              <- max(150,sum(!is.na(Pmat) & Pmat == 1)/2)
	var_i           <- rep(NA,NN)
	var_i[1]        <-  varblock(Pmat)
	for (i in 1:NN){
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
		#pd              <- diff(rbind(Pmat, interval))
		# eligible, 1
		#on_bottom <- pd == -1
		# eligible 2 (on left)
		#on_left   <- t(diff(rbind(0,t(Pmat)))) == 1
		# eligible cells
		#from_cells <- on_bottom & on_left & (prob_from > 0)
		from_cells <- Pmat == 1 & (prob_from > 0)
		#row_from  <- which.max(rowSums(from_cells,na.rm=TRUE) * prob_from)
	    probs      <- rowSums(from_cells, na.rm = TRUE) * prob_from
		row_from  <- sample(1:radix, size = 1, prob = rowSums(from_cells, na.rm = TRUE) * prob_from)
		
		# which has the most zero spots in the to-range?
		col_prob <- colSums(Pmat == 0 * prob_to, na.rm=TRUE) * (from_cells[row_from, ] & !is.na(from_cells[row_from, ]))
		col_in    <- sample(1:n,size=1,prob=col_prob)
		# to cells are zeros in the to-range
		#row_to   <- which.max(Pmat[,col_in] == 0 & (prob_to > 0) * prob_to)
		row_to   <- sample(1:radix, size = 1, prob = rowSums(Pmat[,col_in,drop=FALSE] == 0,na.rm=TRUE) * prob_to)
		# make the change
		Pmat[row_from, col_in] <- 0
		Pmat[row_to, col_in]   <- interval
		var_i[i+1]           <- varblock(Pmat)
		if (i > 100){
			if ((var_i[i+1] - var_i[i-50]) == 0){
				break
			}
		}
	}
	Pmat
}
lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=.5,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
LTblock3 <- LTblock_min(lx_05, prev_05,radix=500,interval=.5)
pdf("REVES2019/Figs/Min1.pdf")
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
