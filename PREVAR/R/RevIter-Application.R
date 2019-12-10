# parametrize only on adjusted lambdas because of reverse iteration only runs like that: 
# doesn´t fit all ages at the same time

### prev mirror (more positive more compressed)
prev_sv_simetric = function(x, y, S, lambda=3){ # x vector of ages at death
  xS = pmax(x-(y-S), 0)
  p <- (xS * exp(abs(lambda) * (xS - S)))/S
  p[p > 1 | is.nan(p)] <- 0
  if(lambda<0){
    pi = sort(1 - p[p>0])
    p[p>0] = pi
  }
  p
}
### que lambda obtiene prevalencia en años person
get_lambda = function(lambda, S, x_ini, x_fin, PrevObs_x, lives){
  Prev_x = integrate(f = prev_sv_simetric, lower =  x_ini, upper = x_fin, subdivisions = 20000, 
                     y = x_fin, S = S, lambda = lambda)$value * lives
  (PrevObs_x - Prev_x)^2}

### get prevalence
get_Prev <- function(lx, ages, lambdas, S){
      dx <- c(diff(-lx), lx[length(lx)])
      prevs = matrix(0, length(ages), length(ages))
      for (y in max(ages):min(ages)){
        prevs[y+1,(y+1):max(y+2-S,0)] = 
            sapply(y:max(y-S+1,0), 
                   FUN = function(s) 
                     integrate(f = prev_sv_simetric, 
                               lower =  s, upper = s+1,  
                               y = y+1, S = S, 
                               lambda = lambdas[y+1],
                               subdivisions = 20000)$value * dx[y+1])}
      prev=colSums(prevs)/lx
      return(prev)
}

### Mortality
x <- 50:109
lx <- c(100000L, 99326L, 99281L, 99249L, 99226L, 99208L, 99194L, 99182L, 
        99169L, 99157L, 99145L, 99134L, 99120L, 99105L, 99088L, 99064L, 
        99031L, 98982L, 98917L, 98824L, 98722L, 98612L, 98487L, 98363L, 
        98237L, 98103L, 97975L, 97842L, 97712L, 97578L, 97442L, 97305L, 
        97163L, 97020L, 96874L, 96729L, 96577L, 96420L, 96252L, 96073L, 
        95884L, 95687L, 95474L, 95235L, 94985L, 94700L, 94389L, 94045L, 
        93669L, 93262L, 92816L, 92333L, 91822L, 91257L, 90648L, 90007L, 
        89307L, 88560L, 87755L, 86898L, 86000L, 85053L, 84077L, 83009L, 
        81914L, 80708L, 79442L, 78077L, 76647L, 75107L, 73480L, 71750L, 
        69942L, 67991L, 65944L, 63758L, 61422L, 58966L, 56389L, 53674L, 
        50857L, 47886L, 44828L, 41659L, 38487L, 35157L, 31738L, 28336L, 
        24930L, 21576L, 18406L, 15280L, 12542L, 9995L, 7776L, 5883L, 
        4370L, 3156L, 2213L, 1505L, 991L, 631L, 389L, 231L, 133L, 74L, 
        40L, 21L, 10L, 5L, 2L)
lx <- lx[50:109]/100
dx <- c(lx[1:59]-lx[2:60], lx[60])
length(lx);length(dx)

### life_bins
# get the ages to which specified quantiles refer.
life_bins <- function(lx, x, probs =  seq(0,1,by=.05)){
  lx <- lx / lx[1]
  splinefun(x~lx, method = "monoH.FC")(probs)
}

### Parameters
S = 10
lambda_limit = c(0,100)

################## Estimate

lambda_estimate <- function(Pi_Obs, lx, S = 10, lambda_limit = c(0,100)){
  
  dx <- c(lx[1:(length(lx)-1)]-lx[2:length(lx)], lx[length(lx)])
  prevs = matrix(0, length(0:59), length(0:59))
  lambda = c()
  # loop for lx
  for (y in 59:0){
    # y = 49
    lambda[y+1] <-optimize(get_lambda, 
                           S = S, 
                           lower = lambda_limit[1], 
                           upper = lambda_limit[2], 
                           x_ini = y, 
                           x_fin = y+1,
                           PrevObs_x = Pi_Obs[y+1] * lx[y+1] - sum(prevs[1:nrow(prevs), y+1]), 
                           lives = dx[y+1])$minimum
    
    prevs[y+1,(y+1):max(y+2-S,0)] = sapply(y:max(y-S+1,0), 
                                           FUN = function(s) 
                                             integrate(f = prev_sv_simetric, 
                                                       lower =  s, upper = s+1,  
                                                       y = y+1, S = S, 
                                                       lambda = lambda[y+1],
                                                       subdivisions = 20000)$value * dx[y+1])
  }
  
  # fit lambdas
  lambdas_df <- data.frame(x,lambda)
  # Variance
  variance <- var(rowSums(prevs))
  return(list(lambdas_df, variance))
}

#################################Application##############################

# fit lambdas adjusted and linearly interpolate

x = 50:109
Pi_Obs1 = 1/(1+exp(-.05*(x-140)))
Pi_Obs4 = 1/(1+.002*exp(-.206*(x-140))) + .01
plot(Pi_Obs1, t="l");lines(Pi_Obs4,lty=2); abline(v=50,lty=3)
Prevs_sim <- sapply(seq(0,1,.05), function(s) Pi_Obs4 + (Pi_Obs1-Pi_Obs4) * s) 
plot(Pi_Obs1,pch=15); abline(v=50,lty=2); for(i in 1:ncol(Prevs_sim)) lines(Prevs_sim[,i], col=i); points(Pi_Obs4,pch=15)

lambdas_adj <- sapply(1:ncol(Prevs_sim), function(s) lambda_estimate(Prevs_sim[,s], lx)[[1]])
variances <- sapply(1:ncol(Prevs_sim), function(s) lambda_estimate(Prevs_sim[,s], lx)[[2]])

beta_min_lala <- function(pars = c(scale = 1, shape1 = 10, shape2 = 5), lambda_obs){
                          xbeta   <- seq(.4, 1,length = length(lambda_obs))
                          lambdas <- pars["scale"] * dbeta(xbeta, shape1 =  pars["shape1"], shape2 =  pars["shape2"])
                          sum((lambda_obs[1:50] - lambdas[1:50])^2)
}

pars_fit <- sapply(1:ncol(lambdas_adj), 
                   function(s) optim(pars_init, beta_min_lala, 
                                             lambda_obs = lambdas_adj[,s]$lambda)$par)

lambdas_fit <- sapply(1:ncol(lambdas_adj), 
                   function(s) pars_fit["scale",s] * dbeta(seq(.4, 1, length.out = 60), 
                              shape1 =  pars_fit["shape1",s], 
                              shape2 =  pars_fit["shape2",s]))

getMSE <- function(pars_fit, Prevs_obs = Prevs_sim, lx, S = 10, ages_fit=0:49){
  
  lambdas = sapply(1:ncol(Prevs_obs), function(s) pars_fit["scale",s] * 
                     dbeta(seq(.4, 1, length.out = 60), 
                           shape1 =  pars_fit["shape1",s], 
                           shape2 =  pars_fit["shape2",s]))
  Prev_estim = sapply(1:ncol(Prevs_obs),
                      function(s) get_Prev(lx, ages = 0:59, lambdas = lambdas[,s], S))
  MSE = sapply(1:ncol(Prevs_obs),
               function(s) sum((Prev_estim[ages_fit,s]-Prevs_obs[ages_fit,s])^2)/
                 length(ages_fit))
  return(list(MSE, Prev_estim))
}

Prev_fit <- getMSE(pars_fit, Prevs_obs = Prevs_sim, lx, S = 10, ages_fit=0:49)



##### graphs
cols <- colorRampPalette(c("red", "blue"))(ncol(Prevs_sim))
par(mfrow=c(1,2))
# prevs to fit
plot(Pi_Obs1,pch=15, cex=.5, main="prevalence profiles"); abline(v=50, lty=2)
for(i in 1:ncol(Prevs_sim)) lines(Prevs_sim[,i], col=cols[i]); points(Pi_Obs4,pch=15,cex=.5)
# lambdas & beta fit
plot(NA,xlim=c(0,60),ylim=c(0,7), main="Lambdas adjusted and Beta fit", ylab="", xlab="age"); abline(v=50, lty=2)
for(i in 1:ncol(lambdas_adj)) points(lambdas_adj[,i]$lambda, col=cols[i], cex=.5, pch=1)
for(i in 1:ncol(lambdas_fit)) lines(lambdas_fit[,i], col=cols[i], lwd=.1)
# variances
plot(NA,xlim=c(0,25),ylim=c(0,100), main="LifeSpan variances", ylab="", xlab="age")
for(i in 1:ncol(lambdas_adj)) points(i, variances[i], col=cols[i], cex=.5, pch=15)
# prevs adj with beta
plot(Pi_Obs1,pch=15, cex=.5, col="white",main="Prevalence Fits",xlab="Age",ylab="Prev")
abline(v=50, lty=2)
for(i in c(1,5,10,15,20)) lines(Prevs_sim[,i], col=cols[i])
for(i in c(1,5,10,15,20)) lines(Prev_fit[[2]][,i], col=cols[i], lty=2)
# errros
plot(1,Prev_fit[[1]][1], xlim=c(0,15), col="white", main="Errors")
for(i in 1:ncol(pars_fit)) points(i,Prev_fit[[1]][i], col=cols[i], pch=15)
# mode (following beta distribution)
modes = (pars_fit[2,]-1)/(pars_fit[2,]+pars_fit[3,]-2)
general_prev = sapply(1:ncol(Prevs_sim), function(i) sum(Prevs_sim[,i]*lx)/sum(lx))
plot(1, modes[1], main="Modes & General Prevalence", ylim=c(.7,1.4), xlim=c(0,.05), 
     xlab="General Prevalence", ylab="Mode")
for(i in 1:ncol(lambdas_fit)) points(general_prev[i],modes[i],col = cols[i], pch=15) 
# shapes
plot(NA,xlim=c(0,15),ylim=c(0,5), main="Shapes", ylab="2", xlab="1")
for(i in 1:ncol(pars_fit)) points(pars_fit[2,i], pars_fit[3,i], col=cols[i], pch=15)
lines(pars_fit[2,], pars_fit[3,])
par(mfrow=c(1,1))


##### Conclusiones
  # ttd prevalences can be fitted with lambda model
  # lambdas by age can be modeled with a beta distribution
  # more Prev compression --> less Beta compression, less mode, more symmetric
  # you tell me Prevalence --> I tell you beta shapes --> variance healthy expectancy 
