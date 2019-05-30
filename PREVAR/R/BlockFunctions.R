
# Author: tim
###############################################################################


makeLTblock <- function(lx, prev, radix, interval){
	n       <- length(lx)
	stopifnot(length(prev) == n)
	lx      <- round(lx / lx[1] * radix)
	Nsick   <- round(lx * prev)
	LTBlock <- matrix(0, ncol = length(lx), nrow = radix)
	for (i in 1:n){
		if (Nsick[i] > 0){
			ind <- sample(1:lx[i], size = Nsick[i], replace = FALSE)
			LTBlock[ind,i] <- interval
		}
	}
	LTBlock
}

# variance simulation
vari <- function(lx, prev, radix = 1e5, interval=1){
	LTBlock <- makeLTblock(
			lx = lx, 
			prev = prev, 
			radix = radix, 
			interval = interval)
	Di      <- rowSums(LTBlock)
	DLE     <- mean(Di) # same every time!
	mean((Di - DLE)^2)
}

#plotLTblock <- function(LTblock,age,flip=FALSE,add = FALSE,...){
#	ypos <- row(LTblock) - 1
#	if (flip){
#		ypos <- nrow(LTblock) - ypos
#	}
#	xpos <- col(LTblock)
#	
#	if (!add){
#		xlim <- c(0,max(xpos))
#		ylim <- c(0,max(ypos))
#		plot(NULL,type="n",ann=FALSE,axes=FALSE,xlim=xlim,ylim=ylim)
#	}
#	
#}
# lx <- lx_5; prev<- prev_5; radix <- 100

blocke0 <- function(lx,prev, interval=1){
	dx   <- lx2dx(lx )
	prev <- prev * interval
	pix  <- cumsum(prev) - prev / 2
	Di   <- dx * pix
	sum(Di) / sum(dx)
}

blockv0 <- function(lx, prev, interval = 1){
	dx   <- lx2dx(lx )
	prev <- prev * interval
	pix  <- cumsum(prev) - prev / 2
	Di   <- dx * pix
	e0   <- sum(Di) / sum(dx)
	sum(((pix - e0)^2) * dx) / sum(dx)
}
e0block <- function(LTblock){
	Di      <- rowSums(LTblock)
	mean(Di)
}
varblock <- function(LTblock){
	Di      <- rowSums(LTblock,na.rm=TRUE)
	DLE     <- mean(Di) # same every time!
	mean((Di - DLE)^2)
}
plot_traj <- function(traj,y,interval = 5,...){
	trle <- rle(traj)
	drop <- trle$values == 0 | is.na(trle$values)
	if (all(drop)){
		return(NULL)
	}
	n    <- length(trle$lengths)
	left <- c(0,cumsum(trle$lengths*interval))[-c(n+1)]
	left <- left[!drop]
	w    <- trle$lengths[!drop] * interval
	rect(left,y,left+w,y+1,border = NA, ...)
}
lx_block_line <- function(lx,a,radix=100,...){
	n   <- length(a)
	lx  <- round(lx / lx[1] * radix)
	a2  <- c(a[1],rep(a[-1],each=2),a[n])
	lx2 <- c(rep(lx,each=2))
	lines(a2,lx2,...)
}
plot_LTblock <- function(LTblock, a, lx, interval = 1, radix = nrow(LTblock), add = TRUE,...){
	lx <- round(lx / lx[1] * radix)
	if (!add){
		plot(a, lx, type = 's', las = 1)
	}
	
	for (i in 1:nrow(LTblock)){
		plot_traj(LTblock[i,],y=i-1,interval = interval,...)
	}
}

# absurd cases:
vblock2 <- function(lx, prev, radix = 1e5, interval = 1, plot = FALSE){
	lx    <- round(lx / lx[1] * radix)
	
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		Pmat[(radix-Nsick[i]):radix,i] <- interval
	}
	if (plot){
		a     <- cumsum(rep(interval,length(lx))) - interval
		plot_LTblock(Pmat[radix:1,],lx=lx,a=a,radix=radix,interval=interval,add=FALSE)
	}
	Di <- rowSums(Pmat)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}

#vblock2(lx_5,prev_5,radix=100,interval=5,plot=TRUE)

vblock3 <- function(lx, prev, radix = 1e5, interval = 1, plot = FALSE){
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
	if (plot){
		a     <- cumsum(rep(interval,n)) - interval
		plot_LTblock(Pmat[radix:1,],lx=lx,a=a,radix=radix,interval=interval,add=FALSE)
	}
	Di <- rowSums(Pmat,na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}
# vblock3(lx_5,prev_5,radix=100,interval=5,plot=TRUE)

# Takes average of fixed rewards and top-down prev
vblock4 <- function(lx, prev, radix = 1e5, interval = 1){
	# take arithmetic average of uniform prevalence
	# and drop-from-lx prevalence. Still constrained.
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		Pmat[(ind+1):(ind+Nsick[i]),i] <- interval
	}
	
	Pmat2 <- matrix(prev * interval,nrow=radix,ncol=n,byrow=TRUE)
	Pmat3 <- (Pmat2 + Pmat) / 2
	
	
	Di <- rowSums(Pmat3, na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}


vblock6 <- function(lx,prev,radix=1e5){
	lx    <- round(lx / lx[1] * radix)
	n     <- length(lx)
	Nsick <- round(lx*prev)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		ind  <- radix-lx[i] 
		if (ind > 0){
			Pmat[1:ind,i] <- NA
		}
		if (ind < radix){
			xnew <- (ind+1):radix
			ref  <- min(prev[i],1-prev[i])
			pl   <- seq(from = prev[i]+ref, to = prev[i]-ref,length=lx[i])
			plm  <- mean(pl)
			pl   <- pl * prev[i]/plm
			Pmat[(ind+1):radix,i] <- pl
		}
	}
	Di <- rowSums(Pmat, na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}

# fixed type- constrained in age - linear dropping from max
# at surv toward min at longest lives. Could be better with
# a min prev maybe? Or an exp drop from surv, such that:
# \int exp fun / lx = prev. Similar to previous tries, 
# but in opposite direction.
vblock7 <- function(lx,prev,radix=1e5, interval=1){
	lx    <- round(lx / lx[1] * radix)
	n     <- length(lx)
	Nsick <- round(lx*prev)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		xnew <- (ind+1):radix
		ref  <- min(prev[i],1-prev[i])
		pl   <- seq(from = prev[i]-ref, to = prev[i]+ref,length=lx[i]) * interval
		plm  <- mean(pl)
		pl   <- pl * prev[i]/plm
		Pmat[(ind+1):radix,i] <- pl
	}
	Di <- rowSums(Pmat, na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}

#prev_sv_block_column = function(q, y, lambda=.3){ # x vector of ages at death
#	p        <- (q * exp(lambda * (q - y)))/y
#	p[p > 1 | is.nan(p)] <- 0
#	p
#}
#mean(prev_sv_block_column(0:100,y=100,lambda=.1))
#plot(0:100,prev_sv_block_column(0:100,y=100,lambda=.3))

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
	NN              <- max(150,sum(!is.na(Pmat) & Pmat == interval)/2)
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
		#prob_to[1:from_i]  <- 0
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
		from_cells <- Pmat == interval & (prob_from > 0)
		#row_from  <- which.max(rowSums(from_cells,na.rm=TRUE) * prob_from)
		probs      <- rowSums(from_cells, na.rm = TRUE) * prob_from
		if (all(probs == 0)){
			break
		}
		row_from   <- sample(1:radix, size = 1, prob = probs)
		
		# which has the most zero spots in the to-range?
		
		col_prob  <- colSums(Pmat == 0 * prob_to, na.rm=TRUE) * (from_cells[row_from, ] & !is.na(from_cells[row_from, ]))
		if (all(col_prob==0)){
			break
		}
		col_in    <- sample(1:n,size=1,prob=col_prob)
		# to cells are zeros in the to-range
		#row_to   <- which.max(Pmat[,col_in] == 0 & (prob_to > 0) * prob_to)
		to_probs <- rowSums(Pmat[,col_in,drop=FALSE] == 0,na.rm=TRUE) * prob_to
		if (all(to_probs==0)){
			break
		}
		row_to   <- sample(1:radix, size = 1, prob = to_probs)
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
	cat("\n",i)
	Pmat
}
