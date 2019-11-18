# parametrize only on adjusted lambdas because of reverse iteration only runs like that: 
# doesn´t fit all ages at the same time
<<<<<<< HEAD

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

=======

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

>>>>>>> 56b320f3fde71942cdea2be32e45ddcafbffdf60
### Mortality
x <- 50:99
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
lx <- lx[50:99]/100
dx <- c(lx[1:49]-lx[2:50], lx[50])
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

lambda_estimate <- function(Pi_Obs, S = 10, lambda_limit = c(0,100)){
  
  prevs = matrix(0, length(0:49), length(0:49))
  lambda = c()
  # loop for lx
  for (y in 49:0){
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
                                           FUN = function(s) integrate(f = prev_sv_simetric, lower =  s, upper = s+1,  
                                                                       y = y+1, S = S, 
                                                                       lambda = lambda[y+1],
                                                                       subdivisions = 20000)$value * dx[y+1])
  }
  
  # fit lambdas
  lambdas_df <- data.frame(x,lambda)
  
  # exponential regression
  model_exp <- nls(lambda ~ alpha * exp(beta * x) , data = lambdas_df, 
                   start = list(alpha = 0.01, beta = 0.05))
  lambdas_hat_exp <- predict(model_exp, list(x = lambdas_df$x))
  
  # polinomial reg
  model_poli <- lm(data = lambdas_df, lambda ~ poly(x,3))
  lambdas_hat_poli <- predict(model_poli, list(x = lambdas_df$x))
  
  # loess
  model_loess <- loess(data = lambdas_df, lambda ~ x, span = .1)
  lambdas_hat_loess <- predict(model_loess, newdata = lambdas_df$x)
  
  # graf
  # plot(lambdas_df)
  # lines(x, lambdas_hat_exp, col = 'skyblue', lwd = 3)
  # lines(x, lambdas_hat_poli, col = 'pink', lwd = 3)
  # lines(x, lambdas_hat_loess, col = 'violet', lwd = 3)
  # legend("topleft", c("exp","cubic", "loess"), lty=2, col= c('skyblue','pink', "violet"), bty = "")
  
  # Prev
  prevs_exp = prevs_loess = matrix(0, length(0:49), length(0:49))
  for (y in 49:0){
    prevs_exp[y+1,(y+1):max(y+2-S,0)] = 
      sapply(y:max(y-S+1,0), FUN = function(s) integrate(f = prev_sv_simetric, 
                                                         lower =  s, upper = s+1,  
                                                         y = y+1, S = S, 
                                                         lambda = lambdas_hat_exp[y+1],
                                                         subdivisions = 20000)$value * dx[y+1])
    prevs_loess[y+1,(y+1):max(y+2-S,0)] = 
      sapply(y:max(y-S+1,0), FUN = function(s) integrate(f = prev_sv_simetric, 
                                                         lower =  s, upper = s+1,  
                                                         y = y+1, S = S, 
                                                         lambda = lambdas_hat_loess[y+1],
                                                         subdivisions = 20000)$value * dx[y+1])
  }
  
  #### results
  par(mfrow=c(1,2))
  plot(x, lambda, main="Lambdas", col = 2, t="p", pch=16, cex=.5, ylim=c(0,5))
  lines(x, lambdas_hat_exp, col=4, lwd=2)
  lines(x, lambdas_hat_loess, col=6, lwd=2)
  legend("topleft", c("Adjusted","Exp","Loess"), lty=1, col= c(2,4,6), bty = "n")
  
  plot(x, Pi_Obs, t = "p", main="Estimate vs Observed",  ylim=c(0,.5),
       xlab="age", ylab= "YP Prev")
  lines(x, round(colSums(prevs)/lx,2), t="l", col=2, lwd=2)
  lines(x, round(colSums(prevs_exp)/lx,2), t="l", col=4, lwd=2)
  lines(x, round(colSums(prevs_loess)/lx,2), t="l", col=6, lwd=2)
  legend("topleft", c("Obs","Adjusted","Exp","Loess"), 
         pch=c(1,NA,NA,NA), lty=c(NA,1,1,1), col= c(1,2,4,6), bty = "n")
  par(mfrow=c(1,1))
  
  return(list(lambdas_df))
<<<<<<< HEAD

}

################### THREE escenarios Prevalence (based in ...)

Pi_Obs1 = 1/(1+exp(-.05*(x-120)))
Pi_Obs3 = 1/(1+.012*exp(-.26*(x-120))) + .025
Pi_Obs2 = (Pi_Obs1 + Pi_Obs3)/2
Pi_Obs4 = 1/(1+.002*exp(-.355*(x-120))) + .02
Pi_Obs5 = (Pi_Obs1+Pi_Obs2)/2


plot(x, lx, t="s", main="lx")
plot(x, Pi_Obs1, t="s", col=2, ylim=c(0,.3), ylab = "Prevalence")
lines(x, Pi_Obs2, t="s", col=4)
lines(x, Pi_Obs3, t="s", col=3)
lines(x, Pi_Obs4, t="s", col=6)
lines(x, Pi_Obs5, t="s", col=1)
title("Time to death Prevalence scenarios")
legend("topleft", c("Prev1","Prev2","Prev3","Prev4"), 
       lty=1, col= c(2,4,3,6), bty = "n")

l1 = lambda_estimate(Pi_Obs1)
l2 = lambda_estimate(Pi_Obs2)
l3 = lambda_estimate(Pi_Obs3)
l4 = lambda_estimate(Pi_Obs4)
l5 = lambda_estimate(Pi_Obs5)
=======

}

################### THREE escenarios Prevalence (based in ...)

Pi_Obs1 = 1/(1+exp(-.05*(x-120)))
Pi_Obs3 = 1/(1+.012*exp(-.26*(x-120))) + .025
Pi_Obs2 = (Pi_Obs1 + Pi_Obs3)/2
Pi_Obs4 = 1/(1+.002*exp(-.355*(x-120))) + .02
Pi_Obs5 = (Pi_Obs1+Pi_Obs2)/2


plot(x, lx, t="s", main="lx")
plot(x, Pi_Obs1, t="s", col=2, ylim=c(0,.3), ylab = "Prevalence")
lines(x, Pi_Obs2, t="s", col=4)
lines(x, Pi_Obs3, t="s", col=3)
lines(x, Pi_Obs4, t="s", col=6)
lines(x, Pi_Obs5, t="s", col=1)
title("Time to death Prevalence scenarios")
legend("topleft", c("Prev1","Prev2","Prev3","Prev4"), 
       lty=1, col= c(2,4,3,6), bty = "n")

l1 = lambda_estimate(Pi_Obs1)
l2 = lambda_estimate(Pi_Obs2)
l3 = lambda_estimate(Pi_Obs3)
l4 = lambda_estimate(Pi_Obs4)
l5 = lambda_estimate(Pi_Obs5)

plot(x, l1[[1]]$lambda, col=2, t="o", ylim = c(0,5))
lines(x, l2[[1]]$lambda, col=4, t="o")
lines(x, l3[[1]]$lambda, col=3, t="o")
lines(x, l4[[1]]$lambda, col=6, t="o")
lines(x, l5[[1]]$lambda, col=1, t="o")

legend("topleft", c("Prev1","Prev","Prev3", "Prev4"), lty=1, col= c(2,4,3,6), bty = "n")
title("Lambda shape related with time to death prevalence profile")

lll <- data.frame(x=x, prev=l3[[1]]$lambda) 
model <- glm(data = lll, formula = prev ~ x, family = "Gamma")
predicted_lambda <-  predict(model, newdata = as.data.frame(x))
plot(x,lll);lines(x,predicted_lambda, col=2)
>>>>>>> 56b320f3fde71942cdea2be32e45ddcafbffdf60

plot(x, l1[[1]]$lambda, col=2, t="o", ylim = c(0,5))
lines(x, l2[[1]]$lambda, col=4, t="o")
lines(x, l3[[1]]$lambda, col=3, t="o")
lines(x, l4[[1]]$lambda, col=6, t="o")
lines(x, l5[[1]]$lambda, col=1, t="o")
legend("topleft", c("Prev1","Prev","Prev3", "Prev4"), lty=1, col= c(2,4,3,6), bty = "n")
title("Lambda shape related with time to death prevalence profile")

## comments
# extender a 120
# beta dist + scale deterministic
# given 2 parameters beta
# three parameters to cover all space