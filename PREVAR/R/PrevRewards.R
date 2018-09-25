
# Author: tim
###############################################################################
# Function modified from original received in 2016 IDEM course
#' calculate the mean, variance, sd, cv, and skewness of lifespan or state occupancy using a rewards framework.
#' @description Calculate various lifetable quantities using a prevalence rewards framework as proposed by Caswell and Zarulli (2018). Choose between Bernoulli or fixed rewards. Fixed means uniform rewards within age (decimal), whereas Bernoulli means 0,1 rewards within age.

#' @details There is no matrix inversion here. The Fundamental matrix is made using lifetable methods. If prevalence refers to good health and you want statistics of poor health, then specify \code{prevx} as the complement.
#' @references
#' Caswell, H., & Zarulli, V. (2018). Matrix methods in health demography: a new approach to the stochastic analysis of healthy longevity and DALYs. Population Health Metrics, 16, 8. http://doi.org/10.1186/s12963-018-0165-5
#' 
#' @param qx numeric vector of onditional death probabilities
#' @param prevx numeric vector of prevalence at each age
#' @param closeout logical. \code{TRUE} if you want half a reward on death. \code{FALSE} gives full reward.
#' @param type integer. 1 for total lifespan statistics. Anything else will do stats of prevalent state occupancy. Default state occupancy.
#' @param interval. Size of time steps. 
#' @param rewards. Either \code{"fixed"} or \code{"Bernoulli"}. 
#' @return \code{data.frame} of some extra lifetable columns: \code{ex}, \code{var}, \code{st}, \code{cv}, and \code{skew}.

PrevRewards <- function(
		qx,            # qx vector
		prevx,         # prevalence vector
		closeout = TRUE,  # better handling of rewards in moment of death
		type = 2,      # 
		interval = 1,  # otherwise it's prevalence
		rewards = "fixed"){ # "bernoulli" different
	# just a question of where to put it if the time interval is e.g. .5
	s       <- length(qx) # nr of ages
	s1      <- s + 1
	stopifnot(s == length(prevx))
	
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
	R1[s1, s1]                      <- 0
	# TR: modify to give 1/2 reward
	R2[row(R2) - 1 == col(R2)]      <- m2
	R2[s1, s1]                      <- 0
	
	R3[row(R3) - 1 == col(R3)]      <- m3
	R3[s1, s1]                      <- 0
	# alternatively, assume fixed rewards:
	if (closeout){
		R1[nrow(R1), 1:s]           <- m1 / 2
		R2[nrow(R2), 1:s]           <- m2 / 2
		R3[nrow(R3), 1:s]           <- m3 / 2
	}
	
	#Z   <- cbind(diag(rep(1,s)), rep(0, s)) # truncation matrix
	#e   <- rep(1, s1) #vector of ones (112x1)
	
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
    # rho1 <- t(N) %*% (Z %*% ((t(P * R1)) %*% e))
	# faster
	rho1    <- colSums( N * colSums(P * R1)[-s1])
	# TR: very similar
	# colSums(solve(I - U) * m1)
	# solve(I - U)
#    .rho2 <- NN %*% 
#		(Z %*% 
#			((t(P * R2)) %*% e) + 
#			2 * Z %*% t(P*R1) %*% t(Z) %*% rho1
#			) 
#    .rho2 <- NN %*%
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
    # almost, but misses closeout
	#colSums(N * (px * m2 + colSums(2 * px * m1 * N * m1)))
	# match
	#.rho2 - rho2
	# holy shit this can be even faster, no need for matrices hardly at all.
	#rho2 <- colSums(N * ((px * m2) + 2*px*m1*rho1))
#	   rho3 <-  NN %*% (                
#				   (t(P * R3) %*% e)[-s1,] +
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
	# similar, but closeouts not handled.
	#colSums(N * c(px*m3 + 3 * px * m2 * rho1 + 3 * px * m1 * rho2))
# original Caswell code:
#rho1 = solve(I-t(U))%*%(Z%*%((t(P*R1))%*%e))
#rho2 = solve(I-t(U))%*%(Z%*%((t(P*R2))%*%e)+2*Z%*%(t(P*R1))%*%t(Z)%*%rho1)
#rho3 = solve(I-t(U))%*%(Z%*%((t(P*R3))%*%e)+3*Z%*%(t(P*R2))%*%t(Z)%*%rho1 + 3*Z%*%(t(P*R1))%*%t(Z)%*%rho2);

	# Statistics of Longevity
	.var        <- rho2 - rho1^2
	# closeout?
    if (.var[s] < 0){
		.var[s] <- 0
	}
	
	#
	std         <- sqrt(.var)
	cv          <- std / rho1
	
	# close out var carefully for sake of skewness
	vv          <- .var
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
	data.frame(ex = rho1, var = .var, st = std, cv = cv, skew = skew)
}


