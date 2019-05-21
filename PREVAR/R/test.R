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
lambda_ex = seq(0, 50, .05)
col<-rainbow(length(lambda_ex))
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,10), xlab='x', ylab='prev')
for (i in 1:length(lambda_ex)){
  lines(seq(0, 10, .05), prev_sv(seq(0, 10, .05), 10, lambda_ex[i]), col=col[i])
}

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

optim_funct = function(y, Pi_obs, iter = 10){
  
  # number individuals
  n = length(y) 
  
  # integer ages of exposures
  Ages = 0:trunc(max(y))
  
  # person-year survivors in each integer age (decimal y gives the portion also)
  lives = sapply(Ages, function(x) sum(pmin(pmax(y-x, 0), 1)))
  
  # obj function
  fun_opt2 = function(lambda, Ages, y, lives, Pi_obs){
    Pi_hat_Age = c()
    for (age in Ages){
      # age=69
      Pi_hat_Age_i = sum(sapply(1:n, function(i) integrate(f = prev_sv, lower =  age, upper = age+1,
                                                           y = y[i], lambda = lambda[i])$value))
      Pi_hat_Age = c(Pi_hat_Age, Pi_hat_Age_i)
    }
    sum((Pi_obs - Pi_hat_Age/lives)^2)
  }
  
  # ineq function
  fun_ineq = function(lambda, Ages, y, lives, Pi_obs){
    z = list()
    for (i in 1:(length(lambda)-1)){
      z[[i]] = lambda[i]-lambda[i+1]
    }
    return(unlist(z))
  }
  
  # iter optim function
  require(Rsolnp)
  Prev_hat = matrix(0, nrow=1, ncol=length(Ages))
  lambda_hats = matrix(0, nrow=1, ncol=n)
  
  for (s in 1:iter){
    
    # initial values
    lambda0 = sample(seq(1,10,.05),n)
    
    # optim
    ineq_bounds = c(0, 10)
    param_bounds = c(0, 50)
    optim = solnp(lambda0, fun = fun_opt2, ineqfun = fun_ineq,
                  ineqLB = rep(ineq_bounds[1], n-1), ineqUB =  rep(ineq_bounds[2], n-1), 
                  LB = rep(param_bounds[1], n), UB = rep(param_bounds[2], n),
                  control = list(delta = .01, rho = 1.5, outer.iter = 50), 
                  Ages = Ages, y = y, lives = lives, Pi_obs = Pi_obs$Pi)
    lambda_hat = matrix(optim$pars, nrow=1)
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
lives_lx <- life_bins(lx, 0:110, probs = seq(1,0,by=-.1))[-1]
plot(lives_lx, seq(1,.01,by=-.1), type = 's')

# give some prevalence
Ages = 0:trunc(max(lives_lx))
Pi_obs_exp = data.frame(x = Ages, Pi = .005*exp(Ages*.047)) # Exponential
Pi_obs_log = data.frame(x = Ages, Pi = .9/(1+exp(-.1*(Ages-80)))) # Logistic
plot(Pi_obs_exp, t='l')
lines(Pi_obs_log, col=2)
Pi_obs = Pi_obs_log # selection

# get 'fit' compression
output_optim = optim_funct(y = lives_lx, Pi_obs, iter = 10) # each iter a starint feasible point
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
error.measures
plot(error.measures$Ages, error.measures$Prev_obs, ylab='Prev', xlab='Ages', t='l')
points(error.measures$Ages, error.measures$Prev_estim_mean, ylab='Prev', xlab='Ages', t='o', col = 2)
legend('topleft', bty='n', legend = c('Obs', 'Est (mean)'), col=1:2, lty=1)
plot(error.measures$Prev_obs, error.measures$Prev_estim, ylab='Estimated', xlab='Observated')
abline(a=0,b=1)

