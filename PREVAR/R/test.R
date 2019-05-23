############################ OPtim with constraints
# Ages >= lambdas: mas ecuaciones q incognitas

# Prev function
# x decimal evaluated age
# y decimal death age
# lambda is compression parameter | higher = more compressed
prev_sv = function(x, y, lambda=.3){ # x vector of ages at death
	p = ifelse(is.infinite(x*exp(x*lambda)), 0, x*exp(x*lambda))/(y*exp(y*lambda))
	p = ifelse(p>1,0,p)
	p
}


# show how flexible lambda is
#lambda_ex = seq(0, 50, .05)
#col<-rainbow(length(lambda_ex))
#plot(0,0, col=0, ylim = c(0,1), xlim=c(0,10), xlab='x', ylab='prev')
#for (i in 1:length(lambda_ex)){
#  lines(seq(0, 10, .05), prev_sv(seq(0, 10, .05), 10, lambda_ex[i]), col=col[i])
#}

# show how flexible lambda is
#lambda_ex = seq(0, 50, .05)
#col<-rainbow(length(lambda_ex))
#plot(0,0, col=0, ylim = c(0,1), xlim=c(0,10), xlab='x', ylab='prev')
#for (i in 1:length(lambda_ex)){
#  lines(seq(0, 10, .05), prev_sv(seq(0, 10, .05), 10, lambda_ex[i]), col=col[i])
#}

# Example: time spent in disability
lambda_x <- 15
y = 3.2
x <- seq(0, y, by=.05)
ts_dis_x = y - splinefun(x~prev_sv(x, y, lambda_x), method = "monoH.FC")(seq(1,0,by=-.01))
plot(x, prev_sv(x, y, lambda_x), col = "red", type = 'l')
lines(ts_dis_x, seq(1,0,by=-.01), col = "blue", lty = 2)
legend(y/2,.8, legend = c('prevalence', 'time spent disability'), lty=c(1,2), col=c(2,4), bty='n')
mean_x = mean(ts_dis_x)
var_x = sum((ts_dis_x - mean_x) ^ 2)
cv_x = var_x^.5/mean_x # tailed distribution

###-------------------------------------------- General Optim function-------------------------------------------
# obj function
fun_opt2 = function(lambda, Ages, y, lives, Pi_obs, decreasing = TRUE){
  n = length(y)
  Pi_hat_Age = c()
  for (age in Ages){
    # age=69
    Pi_hat_Age_i = sum(sapply(1:n, function(i) integrate(f = prev_sv, lower =  age, upper = age+1,
                                                         y = y[i], lambda = lambda[i])$value))
    Pi_hat_Age = c(Pi_hat_Age, Pi_hat_Age_i)
  }
  sum((Pi_obs - Pi_hat_Age/lives)^2)

}

## ineq function
#fun_ineq = function(lambda, Ages, y, lives, Pi_obs){
#	z = list()
#	for (i in 1:(length(lambda)-1)){
#		z[[i]] = lambda[i]-lambda[i+1]
#	}
#	return(unlist(z))
#}

fun_ineq = function(lambda, Ages, y, lives, Pi_obs, decreasing = TRUE){
    if (decreasing){
		z <- diff(lambda)
	} else {
		z <- -diff(lambda)
	}
	z
}

# super saturated one

optim_funct = function(y, Pi_obs, iter = 10){
  # iter optim function
  require(Rsolnp)
  # number individuals
  n = length(y) 
  
  # integer ages of exposures
  Ages = 0:trunc(max(y))
  
  # person-year survivors in each integer age (decimal y gives the portion also)
  lives = sapply(Ages, function(x) sum(pmin(pmax(y-x, 0), 1)))
  
  Prev_hat = matrix(0, nrow=1, ncol=length(Ages))
  lambda_hats = matrix(0, nrow=1, ncol=n)
  
  for (s in 1:iter){
    
    # initial values
    lambda0 = sample(seq(1,10,.05),n)
    
    # optim
    ineq_bounds = c(0, .5)
    param_bounds = c(0, 50)
    optim.s = solnp(lambda0,    # starting pars
                    fun = fun_opt2,     # minimize this
                    ineqfun = fun_ineq, # relationship between pars
                    ineqLB = rep(ineq_bounds[1], n-1), 
                    ineqUB =  rep(ineq_bounds[2], n-1), 
                    LB = rep(param_bounds[1], n), 
                    UB = rep(param_bounds[2], n),
                    control = list(delta = .01, rho = 1.5, outer.iter = 50), 
                    Ages = Ages, 
                    y = y, 
                    lives = lives, 
                    Pi_obs = Pi_obs$Pi)

    lambda_hat = matrix(optim.s$pars, nrow=1)
    lambda_hats = rbind(lambda_hats, lambda_hat)
    
    # save Prev calculated
    Pi_hat_Age = c()
    for (age in Ages){
      Pi_hat_Age_i = sum(sapply(1:n, function(i) integrate(f = prev_sv, lower =  age, upper = age+1,
                                                           y = y[i], lambda = lambda_hat[i])$value))
      Pi_hat_Age = c(Pi_hat_Age, Pi_hat_Age_i)
    }
    Prev_hat =  rbind(Prev_hat, Pi_hat_Age/lives)
  }
  Prev_hat = Prev_hat[-1,]
  lambda_hats = lambda_hats[-1,]
  return(list(Prev_hat, lambda_hats))
}


###-----------------------------------toy data--------------------------------------------------
# invent toy data
y = 1:20

# Invent a prevalence data

Pi_obs = data.frame(x = 0:(max(y)-1), Pi = .005*exp(0:(max(y)-1)*.25)) # at begining of ages
plot(0:(max(y)-1), Pi_obs$Pi, xlab = 'Age', ylab = 'Prev', t='s', xlim = c(0, max(y)), ylim = c(0,1))

# get results

output_optim = optim_funct(y, Pi_obs, iter = 5)
Prev_hat = output_optim[[1]]
lambda_hats = output_optim[[2]]

# graph prevalences
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,max(y)), xlab='Age', ylab='prev')
points(Pi_obs$x+.5, Pi_obs$Pi, col=4, lty=2, lwd=2, t='o')
x = seq(0, max(y), .05)
for (j in 1:nrow(lambda_hats)){
  for (i in 1:length(lambda_hats[j,])){
    lines(x[x<=y[i]], prev_sv(x[x<=y[i]], y[i], lambda_hats[j,i]),col=i)
  }
}

# error Measure
Pi_obs_matrix = matrix(rep(Pi_obs$Pi, nrow(lambda_hats)), nrow= nrow(lambda_hats), ncol=max(y), byrow=T)
error.Abs_iter = Prev_hat - Pi_obs_matrix
error.Relat_iter = round(error.Abs_iter/Pi_obs_matrix*100,2)
error.measures = data.frame(Ages = Ages,
                            Prev_obs = round(Pi_obs$Pi, 4), 
                            Prev_estim = round(colMeans(Prev_hat), 4),
                            Error_mean = round(colMeans(abs(error.Abs_iter)), 4),
                            Error_RealativePorc_mean = round(colMeans(abs(error.Relat_iter)),4))

# ------------------------------------- apply to survival curve-------------------------------------

# interpolate age at survival
life_bins <- function(lx, x, probs =  seq(0,1,by=.05)){
  lx <- lx / lx[1]
  splinefun(x~lx, method = "monoH.FC")(probs)
}

# HMD data
library(HMDHFDplus)
us = 'ivanwilliams1985@gmail.com'; pw = 'volveroman'
mlt <- readHMDweb("USA","mltper_1x1", us, pw)
lx <- subset(mlt,Year == 2010, "lx")$lx

# get ages when survival is 1, .99, .98, ..., .01, 0
q <- seq(1,0,by=-.01)
lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
plot(lives_lx, q[-1] , type = 's')

# give some prevalence
Ages = 0:trunc(max(lives_lx))
Pi_obs_exp = data.frame(x = Ages, Pi = .005*exp(Ages*.047)) # Exponential
Pi_obs_log = data.frame(x = Ages, Pi = .9/(1+exp(-.1*(Ages-80)))) # Logistic
plot(Pi_obs_exp, t='l')
lines(Pi_obs_log, col=2)
Pi_obs = Pi_obs_log # selection

# get 'fit' compression
output_optim = optim_funct(y = lives_lx, Pi_obs, iter = 1) # each iter a starint feasible point
Prev_hat = output_optim[[1]]
lambda_hats = output_optim[[2]]

# graph prevalences
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,trunc(max(lives_lx))+1), xlab='Age', ylab='prev')
points(Pi_obs$x, Pi_obs$Pi, col=4, lty=2, lwd=2, t='l')
x = seq(0, max(lives_lx), .05)
for (j in 1:nrow(lambda_hats)){
  for (i in 1:length(lambda_hats[j,])){
    lines(x[x<=lives_lx[i]], prev_sv(x[x<=lives_lx[i]], lives_lx[i], lambda_hats[j,i]),col=i)
  }
}

# error Measure - diagnostic
Pi_obs_matrix = matrix(rep(Pi_obs$Pi, nrow(lambda_hats)), nrow= nrow(lambda_hats), byrow=T)
error.Abs_iter = Prev_hat - Pi_obs_matrix
error.Relat_iter = round(error.Abs_iter/Pi_obs_matrix*100,2)
error.measures = data.frame(Ages = Ages,
                            Prev_obs = round(Pi_obs$Pi, 4), 
                            Prev_estim_mean = round(colMeans(Prev_hat), 4),
                            Error_mean = round(colMeans(abs(error.Abs_iter)), 4),
                            Error_RealativePorc_mean = round(colMeans(abs(error.Relat_iter)),4))
plot(error.measures$Ages, error.measures$Prev_obs, ylab='Prev', xlab='Ages', t='l')
points(error.measures$Ages, error.measures$Prev_estim_mean, ylab='Prev', xlab='Ages', t='o', col = 2)
legend('topleft', bty='n', legend = c('Obs', 'Est (mean)'), col=1:2, lty=1)
plot(error.measures$Prev_obs, error.measures$Prev_estim, ylab='Estimated', xlab='Observated')
abline(a=0,b=1)


# --------------------------------- analytic prevs -----------------------------------------------


# a utility to help optimize for t_star, which scales the prevalence area by 
# shifting the prevalence curve (implied by alpha and lambda) left or right out of the window S
get_shift = function(t_star, alpha, lambda, S){ 
	
	# https://www.wolframalpha.com/input/?i=%5Cint_0%5Et%7Bx%2Ft*%5Cexp((x-t)*%5Clambda)%7D
	# Area_hat = (lambda * t_star + exp(lambda * -t_star) - 1) / (lambda^2) / t_star
	
	if (sign(t_star) == 1){
		Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S, lambda = lambda)$value + t_star
		
	} else {
		if (sign(t_star) == -1){
			Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S + t_star, lambda = lambda)$value
		} 
	}
	
	((alpha * (S/2)) - (Area_lambda_S)) ^ 2
}
# get the within-lifespan prevalence function
# where:
# x is an age grid, not necessarily integer
# y is the age at death, not necessarily integer!! (probably coming from quantiles)
# lambda gives the shape
# alpha scales the area
# S gives the prevalence window pegged to end-of life

# lambda = 0, alpha = 1 gives a triangle are equal to S/2
# lambda = 0, alpha = 2 fills the rectangle.
# lambda = 0, alpha =.5 shifts the triangle right (beyond y) until the area is equal to .5 * S/2

# lambda > 0, alpha = 1 gives an exp curve shifted left until the area is equal to S/2
# lambda > 0, alpha = 2 fills the rectangle S
# lambda > 0, alpha < 1 moves the curve right until area equal to S/2 * alpha

# remember, needs to return prev for quantile on arbitrarily fine x grid, where the first values
# are probably a bunch of 0s for ages > S. Make sense?

# if pad ==TRUE, then pad with 0s on the right out until omega, otherwise only
# return prevalence curve up until y. Make sense? Marginal calcs will be easier
# if all lifespan prev vectors are the same length.

prev_lambda_alpha_y <- function(
		x = seq(0,omega,by=delta), 
		y = max(x), 
		lambda = 0, 
		alpha = 1.5, 
		S = 20, 
		omega = 110,
		delta = .05){

  # start with....
  # 0      y-S      y    t*
  # |______|________|______|
  #            S      t*-S

   t_star = optimize(get_shift, interval = c(-S, S), alpha = alpha, lambda = lambda, S = S, tol = 1e-12)$minimum
  
   t_star <- round(t_star * 1/delta) * delta

  prev   <- x * 0
  if (sign(t_star) == 1){
	  prev_x <- prev_sv(x = seq(t_star, S, by = delta), y = S, lambda = lambda) 
	  
	  x_implied <- seq((y - S), (y - t_star), by = delta)
	  
	  prev[x >= (y - S) & x <= (y - t_star)] <- prev_x[x_implied>=0]
	  prev[x > (y - t_star) & x <= y] <- 1
  }
  if (sign(t_star) == -1){
	  prev_x <- prev_sv(x = seq(0, S + t_star, by = delta), y = S, lambda = lambda)
	  prev[x >= (y-S-t_star) & x <= (y)] <- prev_x
  }
  if (sign(t_star) == 0){
	  prev_x <- prev_sv(x = seq(0, S, by = delta), y = S, lambda = lambda)
	  prev[x >= (y-S) & x <= (y)] <- prev_x
  }
  prev[x > y] <- NA
  
  return(list(x = x, y = prev))
}




# show some different combos
plot(NULL, type = "n", xlim = c(0,50), ylim = c(0,1))
lines(prev_lambda_alpha_y(y=50,alpha=1,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=.5,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=1.5,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=.5,lambda=.5,S=20))
# that is, get the marginal prevalence by age x assuming we specify 
# lambda (single value or vector thereof)
# alpha (single value or vector thereof)
# lx (standard lifetable lx, but inside we'll get quantiles from it)
# the use the above prev_lambda_alpha_y() to get the prev function for each lifespan

prev_x_from_lambda_alpha <- function(
		lx, 
		x = seq(0, omega, by = delta), 
		lambda = .5, 
		alpha = .5, 
		S = 20, 
		delta = .05, 
		omega = 110,
		q = seq(1,0,by=-.01)){
  
  
  # get ages when survival is 1, .99, .98, ..., .01, 0
  
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx[lives_lx>omega] <- omega
  # vectorize to set flexible shapes by age
  lambdas = rep(lambda, length.out = length(lives_lx))
  alphas = rep(alpha, length.out = length(lives_lx))
  Ss = rep(S, length.out = length(lives_lx))
  
  # get prevalence trajectory for each lifespan
  Prev_lifespan = matrix(0, ncol=length(x), nrow=length(lives_lx))
  for (i in 1:length(lives_lx)){
    # i = 100
    Prev_lifespan[i,] = prev_lambda_alpha_y(
			x = x, 
			y = lives_lx[i], 
			lambda = lambdas[i], 
			alpha = alphas[i], 
			S = Ss[i],
			delta = delta,
			omega = omega)$y
   
  }
  Prev_x <- colMeans(Prev_lifespan, na.rm = TRUE)
  # get mean
#  lives = sapply(x, function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
#  Prev_x = colSums(Prev_lifespan)/lives
  return(list(Prev_x = Prev_x, Prev_lifespan = Prev_lifespan, x = x))
}
Prev_lifespan <- prev_x_from_lambda_alpha(lx)$Prev_lifespan

# For fixed rewards the variance is *zero* unless alpha, lambda, S vary over age.
rowSums(Prev_lifespan,na.rm=TRUE) * delta
# for binomal variance we'd need to invert the prevalence curve 
# in each lifespan. Make sense? Or at least do a discrete inversion,
# in case a prevalence doesn't make it to 1 by death.
# that is a to-do.

# -------------------------------------------------
# Now make a demonstrative plot for IW:
q        <-  seq(1, 0, by = -.05)
lives_lx  <- round(life_bins(lx, 0:110, probs = q)[-1])
lives_lx[lives_lx>omega] <- omega
xx        <- seq(0, 110, by = .1)
S         <- lives_lx * .1 # 10% of life assume
alphas    <- seq(.1, .5, length = length(lives_lx))
lambdas   <- seq(5, .5, length = length(lives_lx))


PrevMat <- matrix(ncol = 1101,nrow=20)

plot(NULL, type = 'n', xlab = "Age", ylab = "Proportion", axes = FALSE, xlim = c(0,110), ylim=c(0,1))
rect(0,q[-length(q)],lives_lx,q[-1], border = gray(.8))
lines(c(0,rep(lives_lx,each=2),0), rep(q,each=2))

for (i in 1:length(lives_lx)){
	pxi <- prev_lambda_alpha_y(
			x = xx, 
			y = lives_lx[i], 
			lambda = lambdas[i], 
			alpha = alphas[i], 
			S = S[i],
			delta = .1,
			omega = omega)
	xx <- pxi$x
	px <- pxi$y
	# scale down
	px <- px * .05
	PrevMat[i, ] <- px
	rem <- is.na(px) | px == 0
	ends <- range(xx[!rem])
	polygon(
			x=c(ends[1],xx[!rem],ends[2]),
			y=c(q[i+1], px[!rem] + q[i+1], q[i+1]) ,
			col = "black"
	)
}

prev_L <- rowSums(PrevMat, na.rm=TRUE) 
var_L = sum((prev_L - mean(prev_L)) ^ 2) / length(prev_L)
var_L^.5/mean(prev_L)
# -------------------------------------------------

#prev_lifespan_inverse <- function(x, prev_lifespan,y){
#	# remove NAs
#	nas 			<- is.na(prev_lifespan)
#	prev_lifespan 	<- prev_lifespan[!nas]
#	x 				<- x[!nas]
#	prev_lifespan   <- prev_lifespan / max(prev_lifespan)
#	splinefun(x~prev_lifespan)
#}


#Prev_lifespan_to_prev_dist <- function(Prev_lifespan){
#	# this isn't so straightforward.
#	# we need to get the inverse function from each lifespan prevalance.
#	
#	
#	
#	
#	
#}


Prev_lifespan <- prev_x_from_lambda_alpha(lx)$Prev_lifespan

# a function to bin prevalence into single ages
bin_prev_x <- function(prev_x, 
                       xin = seq(0,110,by=.05), 
                       xout = 0:100){
  delta_implied = diff(xin)[1]
  xout_index = trunc(xin)
  prev_xout = matrix(0, nrow=nrow(prev_x), ncol=length(xout))
  for (ls in 1:nrow(prev_x)){
    age_d = sum(!is.na(prev_x[ls,])) * delta_implied
    for (age in xout){
      # have to take account for deaths between integer ages
      prev_xout[ls,age+1] = sum(prev_x[ls,xout_index==age] * delta_implied) / min(1,age_d-age, na.rm = T)
    }
  }
  prev_xout
}

Prev_lifespan_bin = bin_prev_x(prev_x_from_lambda_alpha(lx)$Prev_lifespan, 
                               xin = seq(0,110,by=.05), 
                               xout = 0:100)

plot(seq(0,110,by=.05),prev_x_from_lambda_alpha(lx, 
		x = seq(0,110,by=.05), 
		lambda = 10,
		alpha = .2,
		S = 20,
		delta= .05,
		omega = 110,
		q = seq(1,0,by=-.002))$Prev_x, type = 'l', ylim = c(0,1))

plot(seq(0,110,by=.05),prev_x_from_lambda_alpha(lx, 
				x = seq(0,110,by=.05), 
				lambda = 0,
				alpha = 1,
				S = 20,
				delta= .05,
				omega = 110,
				q = seq(1,0,by=-.002))$Prev_x, type = 'l', ylim = c(0,1))



# should do similar to above but just spit out the 
# variance statistic.
prevar_from_lambda_alpha <- function(lx, x=0:100, lambda = .5, alpha = .2, S = 1){

  # get ages when survival is 1, .99, .98, ..., .01, 0
  require(e1071)
  
  # get prevs
  prev_x_output = prev_x_from_lambda_alpha(lx, x, lambda, alpha, S)
  Prev_lifespan = prev_x_output[[2]]
  lives = prev_x_output[[3]]
  
  #var measures
  Prevar_x = colSums(Prev_lifespan^2)/lives - (colSums(Prev_lifespan)/lives)^2
  Prevar2_x = sapply(x, function(s) var(tail(Prev_lifespan[,s+1], lives[s+1])))
  Kurtosis_x = sapply(x, function(s) kurtosis(tail(Prev_lifespan[,s+1], lives[s+1])))
  
  return(list(Prevar_x = Prevar_x, Prevar2_x = Prevar2_x, Kurtosis_x = Kurtosis_x))
}


out1 = prev_x_from_lambda_alpha(lx, x=0:100, lambda = .8, alpha = .5, S = 3)
plot(out1$Prev_x)
out2 = prevar_from_lambda_alpha(lx, x=0:100, lambda = .8, alpha = .5, S = 3)
plot(out2$Prevar_x)
points(out2$Prevar2_x)
plot(out2$Kurtosis_x); abline(h=0)

prev_x_from_lambda_alpha <- function(){
	
}

#
prevar_from_lambda_alpha <- function(){
	# should do similar to above but just spit out the 
	# variance statistic.
}

# These two functions should suffice to make the presentation.

# later, if we wanted to optimize lambda and alpha, we could parameterize their pattern over age in a simple way,
# say by having them vary with a line or curve, for four total parameters a,b,a',b' and these would have
# demographic interpretations of their own.



