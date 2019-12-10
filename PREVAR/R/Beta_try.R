library(tidyverse)

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

x <- 0:110

# lives <- life_bins(lx, x, seq(0,1,.005))
# x_lives <- seq(0,1,.005)*110
# lives_50more <- lives[x_lives>=50]
# x_lives_50more <- x_lives[x_lives>=50] - 50
# plot(x_lives, lives, t="l"); lines(50+x_lives_50more, lives_50more, col=2)
# legend("bottomleft", "older50", col=2, lty=1, bty="n")

lives_50more <- lx[x>=50]
x_lives_50more <- 50:110 - 50

# get pi from lambdas
getPiEst <- function(lx = lives_50more, 
                     x  = x_lives_50more, 
                     lambdas = runif(length(lx), 0, 2), 
                     S = 10){
                
  # lx=5:1; x=20:25; omega = 25; lambdas = runif(length(lx),0, 2); S=3
  Prev <- matrix(0, length(lx), length(x))
  dif_x <- diff(x)[1]
  for(i in length(lx):1){
    # i = 5
    Prev[i,] =  sapply(x, function(s) integrate(f = prev_sv_simetric, 
                                                lower = s - dif_x, upper = s,  
                                                y = x[i], S = S, 
                                                lambda = lambdas[i],
                                                subdivisions = 20000)$value)
  }
  Prev = colSums(Prev[,-length(x)]) / lx 
  return(Prev)
}
  
# cost function
beta_min <- function(pars = c(scale = 1, shape1 = 50, shape2 = 10), 
                     lx, x, S = 10, piObs){
  
  xbeta   <- seq(.5, 1, length = length(lx))
  
  # lambdas = pars_init["scale"] * dbeta(xbeta, shape1 =  pars_init["shape1"], shape2 =  pars_init["shape2"])
  lambdas <- pars["scale"] * dbeta(xbeta, shape1 =  pars["shape1"], shape2 =  pars["shape2"])
  
  PiEst   <- getPiEst(lx, x, lambdas, S)
  
  sum((PiEst - piObs)^2)
}


# Prev to get
Pi_Obs1 = 1/(1+exp(-.05*(x-120)))
Pi_Obs1_interp <- predict(loess(Pi_Obs1 ~ x), newdata = 50 + x_lives_50more)
plot(x, Pi_Obs1, t="s", col=2); lines(50 + x_lives_50more, Pi_Obs1_interp, lwd=2)

# optim

# initial values
pars_init = c(scale = 1/2, shape1 = 50, shape2 = 10)
plot(dbeta(seq(0,1,.01), pars_init[2], pars_init[3]) * pars_init[1], t="l")

# get it!
pars_opt <- optim(pars_init, beta_min, 
                             lx = lx[between(x,50,99)],
                              x = x[between(x,50,100)], 
                              S = 10, 
                              piObs = Pi_Obs1_interp[1:50])

# accuracy
pars_beta_optim <- pars_opt$par
lambdas_optim <- pars_beta_optim["scale"] * 
                 dbeta(seq(.5,.9999,length.out = 50), 
                       shape1 =  pars_beta_optim["shape1"], 
                       shape2 =  pars_beta_optim["shape2"])
plot(lambdas_optim)
Prev_est <- getPiEst(lx = lx[between(x,50,99)],
                     x = x[between(x,50,100)],
                     lambdas = lambdas_optim,
                     S = 10)

plot(Pi_Obs1_interp, ylim=c(0,.5))
lines(Prev_est) ### not working
