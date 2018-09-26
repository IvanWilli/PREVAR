
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
Res$ex[1] - DLE
Res$var[1];VDLE
# different interpretation?

asymptoteBernoulli <- function(dx,prev,interval=1,closeout=TRUE){
	dx        <- dx / sum(dx)
	n         <- length(dx)
	Ntraj     <- sum(2^(1:n))
	stopifnot(length(prev) == n)
	
	trajs     <- get_traj(c(0, 1), n, closeout = closeout)
	
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
DLE <- .71
VDLE <- asymptoteBernoulli(dx=c(.2,.5,.3), prev=c(.1,.5,.7),interval = 1, closeout=FALSE)
asymptoteBernoulli(dx=c(.2,.5,.3), prev=c(.1,.5,.7),interval = 1, closeout=TRUE)

qx <- lx2qx(c(1,.8,.3))

VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=FALSE, calcs="faster")
VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=FALSE, calcs="orig")

VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=TRUE, calcs="faster")
VarBernoulli(qx,prev=c(.1,.5,.7),interval = 1, closeout=TRUE, calcs="orig")



