### Prevalence (based in ...)
x <- 50:99
Pi_Obs = 1/(1+exp(-.05*(x-120)))
plot(x, Pi_Obs, t="s", col=2, ylim=c(0,1))

### Mortality & Prevalence
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
lines(x, lx[x]/100000, t="s")
barplot(rev(x), rev(lx), horiz = T, col = 0, space = 0, xlim=c(50,100))
lx <- lx[50:99]/100
dx <- c(lx[1:49]-lx[2:50], lx[50])
length(lx);length(dx)



### que lambda obtiene prevalencia en años person
get_lambda = function(lambda, S, x_ini, x_fin, PrevObs_x, lives){
  Prev_x = integrate(f = prev_sv_simetric, lower =  x_ini, upper = x_fin, subdivisions = 20000, 
                     y = x_fin, S = S, lambda = lambda)$value * lives
  (PrevObs_x - Prev_x)^2}

### prev mirror
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

### Estimate
prevs = matrix(0, length(0:49), length(0:49))
lambda = c()
# Parameters
S = 10
lambda_limit = c(0,50)
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


#### results
par(mfrow=c(2,2))
plot(x, round(colSums(prevs),2), t="o", main="Estimate vs Observed", xlab="age", ylab="YP Prev")
lines(x, Pi_Obs * lx, t = "o", col=2)
legend(70,30,c("Estimate", "Obs"), lty=2,col=1:2, bty = "n")
plot(round(colSums(prevs),2), Pi_Obs * lx, main="Prevalence Y-P. qqplot", 
     xlab="Estim", ylab="Obs", asp = 0)
abline(1,1)
plot(x, round(colSums(prevs),2)/Pi_Obs/lx, log="y", main="Prevalence Y-P. Est/Obs", 
     ylab="", xlab="age", asp = 0)
abline(h = 1)
plot(x, lambda, main="lambdas")
abline(h = lambda_limit, col=2)
par(mfrow=c(1,1))


### 
# lambda parametrico
# graficos verticales en quantiles




