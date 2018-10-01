
# Author: tim
###############################################################################
# jump to 10 year probabilities
setwd("/home/tim/git/PREVAR/PREVAR")

source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")

lx        <- qx2lx(qx)

# now group / graduate this one
lx_10     <- exp(granularize(log(lx),interval=10,method = "monoH.FC"))

qx_10     <- lx2qx(lx_10)
dx_10     <- lx2dx(lx_10) # weights

a_10      <- seq(0,110,by=10)
n         <- length(qx_10)

#plot(a_5, dx_5)
#plot(a_5, lx_5)
#plot(a_5, qx_5)
prev_10    <- prev_line(lx_10,.01,.5)
# the number of possible discrete trajectories
(Ntraj    <- sum(2^(1:n)))

get_traj <- function(x, maxn = 4, closeout = FALSE){
	traj <- list()
	for (i in 1:maxn){
		li <- as.matrix(expand.grid( rep(list(x), i)))
		if (closeout){
			li[li[,i] == 1,i] <- .5
		}
		traj <- c(traj,split(li,1:nrow(li)))
	}
	#names(traj) <- lapply(traj,paste,collapse="")
	
	traj
}
# comment out
trajs <- get_traj(c(0, 1), n, closeout = TRUE)
#saveRDS(trajs,file="Data/traj10_110.RDS")

length(trajs) == Ntraj

#trajs <- readRDS("Data/traj5_100.RDS")

# put it in a matrix
Tmat  <- matrix(NA,nrow = Ntraj, ncol = n)
for (i in 1:length(trajs)){
	traji                    <- trajs[[i]]
	Tmat[i, 1:length(traji)] <- traji
}

print(object.size(Tmat), units = "Mb")

P1                <- matrix(prev_10, nrow = Ntraj, ncol = n, byrow = TRUE)
Pmat              <- 1 - P1
ind1              <- !is.na(Tmat) & Tmat > 0
Pmat[ind1]        <- P1[ind1]
Pmat[is.na(Tmat)] <- 1 # use 1 for rowProd

#saveRDS(Pmat, file = "Data/Pmat5_100.RDS")

# conditional on length of life,
# probability of having a particular
# Bernoulli sequence of health states.
prevprod <- apply(Pmat,1,prod)
#hist(log(prevprod))
# probability of length of life x
dx_weight <- rep(dx_10/1e5,times=c(2^(1:n)))

# probability of life trajectory
trajprob <- prevprod * dx_weight

# test
sum(trajprob)

# total life in disability
Di     <- rowSums(Tmat, na.rm = TRUE) * 10

# probability of each Di:
DiDist <- tapply(trajprob, Di, sum)
sum(Di == 0)
Tmat[Di == 0,]
barplot(DiDist,space=0)

# expectancy:
(DLE    <- sum(Di * trajprob)) 
# check ()
sum(as.numeric(names(DiDist)) * DiDist) - DLE

# asymptotic variance: E[(X - E[X])^2]
(VDLE <- sum((as.numeric(names(DiDist)) - DLE) ^ 2 * DiDist))
# E(X^2) - E(X)^2, identical ...
(sum((as.numeric(names(DiDist)) ^ 2) * DiDist) - DLE^2)


# same expectancy, different VAR
VarBernoulli(qx_10, prev_10,  interval = 10, calcs="faster")
VarBernoulli(qx_10, prev_10,  interval = 10, calcs="orig")
#Res$ex[1] - DLE
#Res$var[1];VDLE
# different interpretation?

get_Tmat <- function(x = c(0,1), n, closeout = FALSE){
	Ntraj     <- sum(2^(1:n))
	trajs     <- get_traj(x, maxn = n, closeout = closeout)
	
	stopifnot(length(trajs) == Ntraj)
#saveRDS(trajs,file="Data/traj10_110.RDS")
	
	#length(trajs) == Ntraj
	
#trajs <- readRDS("Data/traj5_100.RDS")
	
# put it in a matrix
	Tmat  <- matrix(NA,nrow = Ntraj, ncol = n)
	for (i in 1:Ntraj){
		traji                    <- trajs[[i]]
		Tmat[i, 1:length(traji)] <- traji
	}
	Tmat
}


asymptoteBernoulli <- function(dx,prev,interval=1,closeout=TRUE){
	dx        <- dx / sum(dx)
	n         <- length(dx)
	Ntraj     <- sum(2^(1:n))
	stopifnot(length(prev) == n)
	
	Tmat <- get_Tmat(x=c(0,1),n=n,closeout=closeout)
	
	P1                <- matrix(prev, nrow = Ntraj, ncol = n, byrow = TRUE)
	Pmat              <- 1 - P1
	ind1              <- !is.na(Tmat) & Tmat > 0
	Pmat[ind1]        <- P1[ind1]
	Pmat[is.na(Tmat)] <- 1 # use 1 for rowProd
	

	prevprod <- apply(Pmat,1,prod)
#hist(log(prevprod))
# probability of length of life x
	dx_weight <- rep(dx,times=c(2^(1:n)))
	
# probability of life trajectory
	trajprob <- prevprod * dx_weight
	
    # test
	#sum(trajprob)
	
# total life in disability
	Di     <- rowSums(Tmat, na.rm = TRUE) * interval
	
# probability of each Di:
	DiDist <- tapply(trajprob, Di, sum)
	
# expectancy:
	DLE    <- sum(Di * trajprob)
# check ()
	#sum(as.numeric(names(DiDist)) * DiDist) - DLE
	
# asymptotic variance: E[(X - E[X])^2]
	sum((as.numeric(names(DiDist)) - DLE) ^ 2 * DiDist)
	
	# identical to:
# sum((Di-DLE)^2*trajprob)
}

#DLE <- .71
#VDLE <- asymptoteBernoulli(dx=c(.2,.5,.3), prev=c(.1,.5,.7),interval = 1, closeout=FALSE)
#asymptoteBernoulli(dx=c(.2,.5,.3), prev=c(.1,.5,.7),interval = 1, closeout=TRUE)
#
#qx <- lx2qx(c(1,.8,.3))
#
#VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=FALSE, calcs="faster")
#VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=FALSE, calcs="orig")
#
#VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=TRUE, calcs="faster")
#VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=TRUE, calcs="orig")
#
DLE <- .71
prev=c(.1,.5,.7)
VarFixed(qx=lx2qx(c(1,.8,.3)),prev=c(.1,.5,.7),interval = 1, closeout=FALSE, calcs="faster")
sum((cumsum(c(.1,.5,.7)) - DLE)^2 * lx2dx(c(1,.8,.3),1))

get_Tmat(n=3,closeout=FALSE)

plot_traj <- function(traj,y,h=1,interval = 1,...){
	trle <- rle(traj)
	drop <- trle$values == 0
	if (all(drop)){
		return(NULL)
	}
	n    <- length(trle$lengths)
	left <- c(0,cumsum(trle$lengths*interval))[-c(n+1)]
	left <- left[!drop]
	w    <- trle$lengths[!drop] * interval
	rect(left,y,left+w,y+h,border = NA, ...)
}

plotTmat <- function(Tmat, flip = TRUE, interval = 1, radix = nrow(Tmat)){
	N <- nrow(Tmat)
	if (flip){
		Tmat <- Tmat[N:1, ]
	}
	
	x <- (col(Tmat) - 1) * interval
	y <- row(Tmat) - 1
	
	rect(x,y,x+interval,y+1,
			col=ifelse(Tmat == 1,"black",NA),
			border=ifelse(Tmat == 0,"black",NA))
	text(x+.5,y+.5,Tmat,col=ifelse(Tmat==1,"white","black"))
}


Tmat2Pmat <- function(Tmat,prev){
	n                 <- ncol(Tmat)
	stopifnot(n==length(prev))
	Ntraj             <- sum(2^(1:n))
	P1                <- matrix(prev, nrow = Ntraj, ncol = n, byrow = TRUE)
	Pmat              <- 1 - P1
	ind1              <- !is.na(Tmat) & Tmat > 0
	Pmat[ind1]        <- P1[ind1]
	Pmat[is.na(Tmat)] <- NA # use 1 for rowProd
	Pmat
}

plotPmat <- function(Pmat, Tmat, flip = TRUE, interval = 1, radix = nrow(Tmat)){
	N <- nrow(Tmat)
	if (flip){
		Tmat <- Tmat[N:1, ]
		Pmat <- Pmat[N:1, ]
	}
	
	x <- (col(Tmat) - 1) * interval
	y <- row(Tmat) - 1
	
	rect(x,y,x+interval,y+1,
			col=ifelse(Tmat == 1,"black",NA),
			border=ifelse(Tmat == 0,"black",NA))
	text(x+.5,y+.5,Pmat,col=ifelse(Tmat==1,"white","black"))
}



plot_Pmat_traj <- function(traj, y, h = 1, interval = 1,
		x = c(0, 1), cols = c(NA, "black"), border = c("black",NA)
){
	trle <- rle(traj)
	drop <- is.na(trle$values)
	
	names(cols)   <- x
	names(border) <- x
	
	colse    <- cols[as.character(trle$values[!drop])]
	borderse  <- border[as.character(trle$values[!drop])]
	if (all(drop)){
		return(NULL)
	}
	n    <- length(trle$lengths)
	left <- c(0,cumsum(trle$lengths*interval))[-c(n+1)]
	left <- left[!drop]
	w    <- trle$lengths[!drop] * interval
	rect(left,y,left+w,y+h,col=colse,border=borderse)
}
plotBernoulliSeq <- function(Pmat, Tmat, interval = 1, flip = TRUE){
	N <- nrow(Tmat)
	if (flip){
		Tmat <- Tmat[N:1, ]
		Pmat <- Pmat[N:1, ]
	}
	hs <- apply(Pmat,1,prod,na.rm=TRUE)
	ys <- cumsum(hs) - hs
	for (i in 1:N){
		plot_Pmat_traj(traj=Tmat[i, ],y=ys[i],h=hs[i],interval = interval,
				x = c(0, 1), cols = c(NA, "black"), border = c("black",NA))
	}
	n <- ncol(Pmat)
	polygon(x=c(0,rep((n:1)*interval,each=2),0),
			y=c(rep(0:n,each=2))
	)
}

plotBernoulliMarkovSeq <- function(Pmat, Tmat, dx, interval = 1, flip = TRUE){
	N <- nrow(Tmat)
	dx_expanded <- rep(dx, times = table(ncol(Tmat) - rowSums(is.na(Tmat))))
	if (flip){
		Tmat <- Tmat[N:1, ]
		Pmat <- Pmat[N:1, ]
		dx   <- rev(dx)
		dx_expanded <- rev(dx_expanded)
	}

	hs <- apply(Pmat,1,prod,na.rm=TRUE) * dx_expanded
	ys <- cumsum(hs) - hs
	for (i in 1:N){
		plot_Pmat_traj(traj=Tmat[i, ],y=ys[i],h=hs[i],interval = interval,
				x = c(0, 1), cols = c(NA, "black"), border = c("black",NA))
	}
	n <- ncol(Pmat)
	polygon(x=c(0,rep((n:1)*interval,each=2),0),
			y=c(rep(c(0,dx),each=2))
	)
}

Tmat <- get_Tmat(n=3,closeout=FALSE)
Pmat <- Tmat2Pmat(Tmat,prev)

# 1) Binary trajectory space
pdf("Figures/BernTraj.pdf",height=7.7,width=2.2)
par(mai=c(.5,.5,.2,.2))
plot(NULL, type = "n", xlim = c(0, 3), ylim = c(0, 14), ann = FALSE, axes = FALSE, asp = 1)
plotTmat(Tmat)
polygon(c(0,3,3,2,2,1,1,0),c(0,0,8,8,12,12,14,14))
axis(1,at=c(0,1,2,3))
text(0,.5:13.5,14:1,pos=2,xpd=TRUE)
dev.off()
# 2) with Bernoulli probabilities
pdf("Figures/BernTrajProbs.pdf",height=7.7,width=2.2)
par(mai=c(.5,.5,.2,.2))
plot(NULL, type = "n", xlim = c(0,3),ylim = c(0,14),ann=FALSE, axes = FALSE,asp=1)
plotPmat(Pmat,Tmat)
polygon(c(0,3,3,2,2,1,1,0),c(0,0,8,8,12,12,14,14))
axis(1,at=c(0,1,2,3))
text(0,.5:13.5,14:1,pos=2,xpd=TRUE)
dev.off()
# 3) with height weighted to probability of trajectory
# ocurring conditional on trajectory length (length of life)
pdf("Figures/BernCondTrajProbs.pdf",height=5.7,width=5.7)
par(mai=c(.5,.5,.2,.2))
plot(NULL, type = "n", xlim =
				c(0,3), ylim = c(0,3), ann = FALSE, axes = FALSE)
plotBernoulliSeq(Pmat,Tmat)
axis(1,at=c(0,1,2,3),cex.axis=1.5)
axis(2,at=c(0,1,2,3),las=1,cex.axis=1.5)
dev.off()
# 4) Bernoulli also Markov weighted
pdf("Figures/BernTrajProbsWeighted.pdf",height=5.7,width=5.7)
par(mai=c(.5,.5,.2,.2))
plot(NULL, type = "n", xlim =
				c(0,3), ylim = c(0,1), ann = FALSE, axes = FALSE)
plotBernoulliMarkovSeq(Pmat,Tmat,dx=c(.2,.5,.3))
axis(1,at=c(0,1,2,3),cex.axis=1.5)
axis(2,at=c(0,1),las=1,cex.axis=1.5)
dev.off()


prevprod <- apply(Pmat,1,prod,na.rm=TRUE)
#hist(log(prevprod))
# probability of length of life x
dx_weight <- rep(c(.2,.5,.3), times=c(2^(1:3)))
ptraj     <- prevprod * dx_weight
Di        <- rowSums(Tmat, na.rm = TRUE)

# probability distribution of total time spent
pdf("Figures/BernDiDist.pdf", height = 5.7, width = 5.7)
par(mai=c(.5,.7,.2,0))
barplot(tapply(ptraj, Di, sum),space=0,las=1,cex.axis=1.5,cex=1.5)
dev.off()