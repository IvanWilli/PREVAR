
# Author: tim
###############################################################################


makeLTblock <- function(lx, prev, radix, interval){
	n       <- length(lx)
	stopifnot(length(prev) == n)
	lx      <- lx / lx[1] * radix
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

# absurd cases:
vblock2 <- function(lx,prev, radix = 1e5){
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		Pmat[(radix-Nsick[i]):radix,i] <- 1
	}
	Di <- rowSums(Pmat)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}

vblock3 <- function(lx,prev, radix = 1e5){
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		Pmat[(ind+1):(ind+Nsick[i]),i] <- 1
	}
	Di <- rowSums(Pmat,na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}

vblock4 <- function(lx, prev, radix = 1e5){
	# take arithmetic average of uniform prevalence
	# and drop-from-lx prevalence. Still constrained.
	lx    <- round(lx / lx[1] * radix)
	Nsick <- round(lx*prev)
	n     <- length(lx)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		Pmat[(ind+1):(ind+Nsick[i]),i] <- 1
	}
	
	Pmat2 <- matrix(prev,nrow=radix,ncol=n,byrow=TRUE)
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
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		xnew <- (ind+1):radix
		ref  <- min(prev[i],1-prev[i])
		pl   <- seq(from = prev[i]+ref, to = prev[i]-ref,length=lx[i])
		plm  <- mean(pl)
		pl   <- pl * prev[i]/plm
		Pmat[(ind+1):radix,i] <- pl
	}
	Di <- rowSums(Pmat, na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}


vblock7 <- function(lx,prev,radix=1e5){
	lx    <- round(lx / lx[1] * radix)
	n     <- length(lx)
	Nsick <- round(lx*prev)
	Pmat  <- matrix(0,nrow=radix,ncol=n)
	for (i in 1:n){
		ind <- radix-lx[i] 
		Pmat[1:ind,i] <- NA
		xnew <- (ind+1):radix
		ref  <- min(prev[i],1-prev[i])
		pl   <- seq(from = prev[i]-ref, to = prev[i]+ref,length=lx[i])
		plm  <- mean(pl)
		pl   <- pl * prev[i]/plm
		Pmat[(ind+1):radix,i] <- pl
	}
	Di <- rowSums(Pmat, na.rm=TRUE)
	e0 <- mean(Di)
	mean((e0 - Di)^2)
}