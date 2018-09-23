# Author: tim
###############################################################################


SullivanMatrixCalc <- function(
        qx,            # qx vector
        prevx,         # prevalence vector
        closeout = TRUE,  # better handling of rewards in moment of death
        type = 1,      # 1 healthy, 2 sick
        interval = 1,  # otherwise it's prevalence
		rewards = "fixed"){ # "bernoulli" different
                       # just a question of where to put it if the time interval is e.g. .5
    s       <- length(qx) # nr of ages
	s1      <- s+1
    stopifnot(length(qx) == length(prevx))


    I       <- diag(s) # identity matrix
    
	# TR: this could be generalized more
	# build transition matrix
    px      <- 1 - qx
    U       <- matrix(0, nrow = s, ncol = s)
    U[row(U)-1 == col(U) ] <- px[-s]
    U[s,s]  <- px[s]

    # define the Markov chain
    P       <- rbind(cbind(U, rep(0, s)), c(1 - colSums(U), 1))
    
	# figure out rewards
	
	# type == 1 for total lifetime
    if(type == 1){
        # reward for living another interval
        m1 <- rep(interval,s)
        m2 <- rep(interval,s)
        m3 <- rep(interval,s)
		
	# otherwise it's prevalence related
    } else {
        # reward 1 (unit of healthy life) with probability 'prevage'
		m1 <- prevx * interval
		if (rewards == "fixed"){
			# fixed means uniform
			m2 <- prevx ^ 2 * interval
			m3 <- prevx ^ 3 * interval
		} else {       
			# Bernoulli outcome
            m2 <- prevx * interval
            m3 <- prevx * interval
		}
    }

    # Define reward matrices with rewards associated with the transitions
    R1 <- matrix(0, nrow = s1, ncol = s1)
    R2 <- R3 <- R1
    # to place entries on the subdiagonal match row to column and then subtract 1 from row
    R1[row(R1) - 1 == col(R1)]      <- m1
    R1[s+1,s+1]                     <- 0
    # TR: modify to give 1/2 reward
    R2[row(R2) - 1 == col(R2)]      <- m2
    R2[s+1,s+1]                     <- 0

    R3[row(R3) - 1 == col(R3)]      <- m3
    R3[s+1,s+1]                     <- 0
    # alternatively, assume fixed rewards:
    if (closeout){
        R1[nrow(R1), 1:s]           <- m1 / 2
        R2[nrow(R2), 1:s]           <- m2 / 2
        R3[nrow(R3), 1:s]           <- m3 / 2
    }

    Z   <- cbind(diag(rep(1,s)), rep(0, s)) # truncation matrix
    e   <- rep(1, s + 1) #vector of ones (112x1)

    #NN   <- solve(I - t(U))
	#N    <- solve(I - U)
    #faster for large s
    N <- matrix(0,s,s)
	for (i in 1:s){
		N[i:s,i] <- qx2lx(qx[i:s],radix=1)
	}
	
    # moments of the accumulated rewards
    # rho1   <- solve(I - t(U)) %*% (Z %*% ((t(P * R1)) %*% e))
    # rho1    <- NN %*% (Z %*% ((t(P * R1)) %*% e))
    # rho1    <- NN %*% ((t(P * R1) %*% e)[-s1,])

	# faster
	rho1    <- colSums( N * colSums(P * R1)[-s1])
    # TR: very similar
    # colSums(solve(I - U) * m1)
    # solve(I - U)

#    rho2 <- NN %*%
#            (
#                (t(P * R2) %*% e)[-s1,] +
#                2 * t(P*R1)[-s1,-s1] %*% rho1
#                )  
#   faster
	rho2 <- colSums(N *
			(
				colSums(P * R2)[-s1] +
				colSums(2 *(P*R1)[-s1,-s1] * c(rho1))
				))

#    rho3 <-  NN %*% (
#                (t(P * R3) %*% e)[-s1,] +
#                3 * t(P * R2)[-s1,-s1] %*% rho1 +
#                3 * t(P * R1)[-s1,-s1] %*% rho2
#                )
	# faster
    rho3 <-  colSums(N *
				(
					colSums(P * R3)[-s1] +
					3 * colSums((P * R2)[-s1,-s1] * c(rho1)) +
					3 * colSums((P * R1)[-s1,-s1] * c(rho2))
					))

    # Statistics of Longevity
    var         <- rho2 - rho1^2
    std         <- sqrt(var)
    cv          <- std / rho1

    # close out var carefully for sake of skewness
    vv          <- var
    pick        <- vv == 0
    vv[pick]    <- 1
    vv          <- vv ^ (-3 / 2)
    vv[pick]    <- 0
    dim(vv)     <- NULL
  
    # TR: left side of equation
   # dd         <- diag(vv, nrow = length(vv))
    #skew      <- dd %*% (rho3 - 3 * rho2 * rho1 + 2 * rho1^3)
    # faster
	skew        <-  vv * (rho3 - 3 * rho2 * rho1 + 2 * rho1^3)
    # save Longevity stats (type 1) in the first 5 columns, and Healthy LE in the next five
    data.frame(ex = rho1, var = var, st = std, cv = cv, skew = skew)
}



lx        <- qx2lx(qx)
# now graduate this one

lx_5      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
lx_10     <- exp(granularize(log(lx),interval=10,method = "monoH.FC"))
lx.5      <- exp(granularize(log(lx),interval=.5,method = "monoH.FC"))
lx.25     <- exp(granularize(log(lx),interval=.25,method = "monoH.FC"))
lx.1      <- exp(granularize(log(lx),interval=.1,method = "monoH.FC"))
lx.025    <- exp(granularize(log(lx),interval=.025,method = "monoH.FC"))
lx.01     <- exp(granularize(log(lx),interval=.01,method = "monoH.FC"))
#lx.5  <- exp(granularize(log(lx),interval=.5,method = "linear"))
#lx.25 <- exp(granularize(log(lx),interval=.25,method = "linear"))
#lx.1  <- exp(granularize(log(lx),interval=.1,method = "linear"))
qx_5      <- lx2qx(lx_5)
qx_10     <- lx2qx(lx_10)
qx.5      <- lx2qx(lx.5)
qx.25     <- lx2qx(lx.25)
qx.1      <- lx2qx(lx.1)
qx.025    <- lx2qx(lx.025)
qx.01     <- lx2qx(lx.01)

prev      <- prev_line(lx,0,.5)
prev.5    <- prev_line(lx.5,0,.5)
prev.25   <- prev_line(lx.25,0,.5)
prev.1    <- prev_line(lx.1,0,.5)
prev.025  <- prev_line(lx.025,0,.5)
prev.01   <- prev_line(lx.01,0,.5)

# how many sick / age
#plot(0:110,qx,log='y',ylim=c(1e-6,1))
#lines(seq(0,110,by=.5),qx.5)
#lines(seq(0,110,by=.25),qx.25)
#lines(seq(0,110,by=.25),qx.25)
#lines(seq(0,110,by=.1),qx.1)
#qx[1]
#1-cumprod(1-qx.1[1:10])
#plot(diff(qx)[2:10])
#diff(qx)

# radiz size, intervals = 1
do.this <- FALSE
if (do.this){
N       <- 1000
Vars1e5 <- rep(0,N)
for (i in 1:N){
	if (i%%100 == 0){
		cat(i,"\n")
	}
	Vars1e5[i] <- vari(lx = lx, prev = prev, radix = 1e5)
}
hist(Vars1e5)

N    <- 1000
Vars1e4 <- rep(0,N)
for (i in 1:N){
	Vars1e4[i] <- vari(lx = lx, prev = prev, radix = 1e4)
}
hist(Vars1e4)

N    <- 1000
Vars1e3 <- rep(0,N)
for (i in 1:N){
	Vars1e3[i] <- vari(lx = lx, prev = prev, radix = 1e3)
}

png("/home/tim/workspace/Other/Figures/VarSim.png")
hist(Vars1e3,breaks = seq(33,39,by=.1),ylim=c(0,600),col="#00FF0080",border=NA,
		xlab = "Variance distributions\n(N = 1000 for each)",
		main = "Simulation results for 3 difference radix sizes")
hist(Vars1e4,add=TRUE,breaks = seq(33,39,by=.1),col="#0000FF80",border=NA)
hist(Vars1e5,add=TRUE,col="#FF000080",breaks = seq(33,39,by=.1),border=NA)
text(36.3,400,paste0("Radix 100k"),pos=2)
text(36.3,160,paste0("Radix 10k"),pos=2)
text(35.5,50,paste0("Radix 1k"),pos=2)
q5 <- quantile(Vars1e5,probs = c(.025,.5,.975))
q4 <- quantile(Vars1e4,probs = c(.025,.5,.975))
q3 <- quantile(Vars1e3,probs = c(.025,.5,.975))
segments(q5[1],400,q5[3],400,lwd=2,col="red")
segments(q4[1],180,q4[3],180,lwd=2,col="blue")
segments(q3[1],70,q3[3],70,lwd=2,col="green")
abline(v=c(q5[2],q4[2],q3[2]),col=c("red","blue","green"))
dev.off()

#abline(v=SullivanMatrixCalc(qx, prev, type = 2, closeout = FALSE)$var[1] )
#abline(v=SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$var[1])
}

# 


# the last one, .01, eats memory, runs over an hour.
#variance_intervals <- c(
#SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$var[1],
#SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval=.5)$var[1],
#SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE, interval=.25)$var[1],
#SullivanMatrixCalc(qx.1, prev.1, type = 2, closeout = TRUE, interval=.1)$var[1],
#SullivanMatrixCalc(qx.025, prev.025, type = 2, closeout = TRUE, interval=.025)$var[1],
#SullivanMatrixCalc(qx.01, prev.01, type = 2, closeout = TRUE, interval=.01)$var[1]
#)
#sd_intervals <- sqrt(variance_intervals)
#dput(sd_intervals)
sd_intervals <- c(5.68979370525243, 5.84416536170791, 5.9214188285196, 5.967797098941, 
		5.99099423181927, 5.99563432234453)
gc()

png("Figures/Convergence.png")
plot(c(1,.5,.25,.1,.025,.01),sd_intervals )
dev.off()
x <- c(1,.5,.25,.1,.025,.01)
mod <- lm(sd_intervals~x)
predict(mod,data.frame(x=0))
summary(mod)

dput(sd_intervals)
variance_fixedrewards <- c(
SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE, rewards = "fixed")$var[1],
SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval=.5, rewards = "fixed")$var[1],
SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE, interval=.25, rewards = "fixed")$var[1],
SullivanMatrixCalc(qx.1, prev.1, type = 2, closeout = TRUE, interval=.1, rewards = "fixed")$var[1],
SullivanMatrixCalc(qx.025, prev.025, type = 2, closeout = TRUE, interval=.025, rewards = "fixed")$var[1])

plot(sd_intervals[1:5],
sqrt(variance_fixedrewards),asp=1)
abline(a=0,b=1)

variance_fixedrewards <- c(
		SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE, rewards = "fixed")$var[1],
		SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval=.5, rewards = "fixed")$var[1],
		SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE, interval=.25, rewards = "fixed")$var[1],
		SullivanMatrixCalc(qx.1, prev.1, type = 2, closeout = TRUE, interval=.1, rewards = "fixed")$var[1],
		SullivanMatrixCalc(qx.025, prev.025, type = 2, closeout = TRUE, interval=.025, rewards = "fixed")$var[1])
variance_hmrewards <- c(
		SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE, rewards = "hm")$var[1],
		SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval=.5, rewards = "hm")$var[1],
		SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE, interval=.25, rewards = "hm")$var[1],
		SullivanMatrixCalc(qx.1, prev.1, type = 2, closeout = TRUE, interval=.1, rewards = "hm")$var[1],
		SullivanMatrixCalc(qx.025, prev.025, type = 2, closeout = TRUE, interval=.025, rewards = "hm")$var[1])

#
#plot(seq(0,110,by=1),SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$var,ylim=c(0,35))
#lines(seq(0,110,by=.5),
#		SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval = .5)$var,
#		col="blue")
#lines(seq(0,110,by=.25),
#		SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE, interval = .25)$var,
#		col = "royalblue")
#lines(seq(0,110,by=.1),
#		SullivanMatrixCalc(qx.1, prev.1, type = 2, closeout = TRUE, interval = .1)$var,
#		col = "magenta")
e0i <- function(lx, prev, radix = 1e5, interval=1){
	n       <- length(lx)
	stopifnot(length(prev) == n)
	lx      <- lx / lx[1] * radix
	Nsick   <- round(lx*prev)
	LTBlock <- matrix(0,ncol = length(lx), nrow = radix)
	for (i in 1:n){
		if (Nsick[i] > 0){
			ind <- sample(1:lx[i],size=Nsick[i],replace=FALSE)
			LTBlock[ind,i] <- interval
		}
	}
	Di      <- rowSums(LTBlock)
	DLE     <- mean(Di) # same every time!
	DLE
}

# ex converges.

e0i(lx.025,prev.025,radix=1e5,interval=.025)
SullivanMatrixCalc(qx.025, prev.025, type = 2, closeout = TRUE,interval=.025)$ex[1]


Pmat <- matrix(prev,nrow=1e5,ncol=length(prev),byrow=TRUE)

lx2dx <- function(lx,radix=1e5){
	lx <- lx / lx[1]
	-diff(radix * c(lx,0))
}
lx2dx(lx)
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
blocke0(lx,prev,1)
blocke0(lx.25,prev.25,.25)

# perfectly uniform out to some degree of prevision
sqrt(blockv0(lx,prev,1))
sqrt(blockv0(lx.25,prev.25,.25))
sqrt(blockv0(lx.01,prev.01,.01))

plot(sqrt(c(blockv0(lx,prev,1),
				blockv0(lx.5,prev.5,.5),
				blockv0(lx.25,prev.25,.25),
				blockv0(lx.1,prev.1,.1),
				blockv0(lx.025,prev.025,.025))),
		sqrt(variance_fixedrewards),asp=1)


N         <- 1000
Vars1e5.1 <- Vars1e5.25 <- Vars1e5.5 <- Vars1e5 <- rep(0,N)
for (i in 1:N){
	Vars1e5[i]    <- vari(lx = lx, prev = prev, radix = 1e5)
	Vars1e5.5[i]  <- vari(lx = lx.5, prev = prev.5, radix = 1e5, interval = .5)
	Vars1e5.25[i] <- vari(lx = lx.25, prev = prev.25, radix = 1e5, interval = .25)
	Vars1e5.1[i]  <- vari(lx = lx.1, prev = prev.1, radix = 1e5, interval = .1)
	if (i %% 10 == 0){
		cat(i,"\n")
	}
}

simulated_sd <- sqrt(c(
mean(Vars1e5),
mean(Vars1e5.5),
mean(Vars1e5.25),
mean(Vars1e5.1)))

plot(sd_intervals[1:4],simulated_sd)

par(mfrow=c(1,2))
plot(c(1,.5,.25,.1),simulated_sd,ylim=c(4.8,6),xlim=c(0,1))
text(c(1,.5,.25,.1),simulated_sd,pos=1,c(1,.5,.25,.1))

plot(c(1,.5,.25,.1,.025,.01),sd_intervals ,ylim=c(4.8,6),xlim=c(0,1))
text(c(1,.5,.25,.1,.025,.01),sd_intervals,pos=1,c(1,.5,.25,.1,.025,.01))
Vars1e3.5 <- rep(0,N)
for (i in 1:N){
	Vars1e4.5[i] <- vari(lx = lx.5, prev = prev.5, radix = 1e4,interval=.5)
}


vS1 <- SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$var[1]
mean(Vars1e4)-vS1
SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$ex[1]
SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval=.5)$ex[1]
# now with interval .5

N         <- 10000
Vars1e4.5 <- rep(0,N)
lx.5      <- granularize(lx,interval=.5)




qx.5      <- lx2qx(lx.5)
prev.5    <- prev_line(lx.5)
for (i in 1:N){
	Vars1e4.5[i] <- vari(lx = lx.5, prev = prev.5, radix = 1e4)
}
vS.5      <- SullivanMatrixCalc(qx.5, prev.5, type = 2, closeout = TRUE, interval = .5)$var[1]
mean(Vars1e4.5)-vS.5

# now with interval .25

N          <- 10000
Vars1e4.25 <- rep(0,N)
Lx.25      <- granularize(Lx,interval=.25)
qx.25      <- lx2qx(Lx.25)
prev.25    <- prev_line(Lx.25)
for (i in 1:N){
	Vars1e4.25[i] <- vari(Lx = Lx.25, prev = prev.25, radix = 1e4)
}
vS.25      <- SullivanMatrixCalc(qx.25, prev.25, type = 2, closeout = TRUE)$var[1]
mean(Vars1e4.25)-vS.25









# absurd case:
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
sqrt(vblock2(lx,prev))
sqrt(vblock3(lx,prev))


barplot(table(Di))

png("/home/tim/git/Spells/Spells/Figures/Top-Down-Bottom-Up.png")
plot(0:110,lx,type='l', xlab = "Age", las = 1, 
		main ="Binary prevalence top-down and bottom up\nlinear increase from 0 at age 0 to .5 at age 110")
polygon(c(0:110,110:0),c(Nsick,rep(0,111)),col="#FF000080",border="red")
polygon(c(0:110,110:0),c(lx-Nsick,rev(lx)),col="#0000FF80",border="blue")
text(60,8e4,"sd = 10.65",cex=2,font=2)
text(60,1e4,"sd = 26.52",cex=2,font=2)
dev.off()
ramp   <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"),space="Lab")
breaks <- seq(0,.5,by=.01)
colors <- as.character(cut(prev,breaks=breaks,labels=ramp(length(breaks)-1)))

png("/home/tim/git/Spells/Spells/Figures/Uniform.png")
plot(0:110,lx,type='l',
		main = "Prevalence Uniform within age\nlinear increase from 0 at age 0 to .5 at age 110",
		las=1,xlab="Age")
for (i in 0:110){
	polygon(c(i,i+1,i+1,i),c(0,0,lx[i+2],lx[i+1]), border= NA, col = colors[i])
}
lines(0:110,lx,lwd=2)
text(40,5e4,"sd= 4.74",cex=2,font=2)
dev.off()
#sqrt(blockv0(lx,prev,1))
# -----------------------------------

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
sqrt(vblock4(lx,prev))
sqrt(SullivanMatrixCalc(qx, prev, type = 2, closeout = TRUE)$var[1])

# and what about top-down linear prevalence?

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
sqrt(vblock6(lx,prev,radix=1e5))
sqrt(vblock7(lx,prev,radix=1e5))


seq(0,.5,length=100)
plot(seq(0,.5,length=100),pexp(q=seq(0,.5,length=100),rate=.01))

lines(0:10,dexp(x=0:10,rate=.2))


# simple to parameterize. we have nr steps to take,
d <- 1/sum(1:lx[10]) # this is the second derivative of the curve
plot(cumsum(cumsum(rep(d,lx[80]))),type='l')

mean(cumsum(cumsum(rep(d,lx[80]))))
