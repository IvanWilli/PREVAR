omega <- 110
delta = .05
q = seq(1,0,by=-.01)
x = seq(0, omega, by = delta)
S <- 20

# select 3 persons at random
lives_lx = c(39.40, 76.20)

# symmetric funct
prev_sv_simetric = function(x, y, lambda=3){ # x vector of ages at death
  lambda_abs = abs(lambda)
  p <- (x * exp(lambda_abs * (x - y)))/y
  if(lambda<0){
    p = sort(1 - p)
  }
  p[p > 1 | is.nan(p)] <- 0
  p
}
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,100), xlab='x', ylab='prev')
for (lambdas in seq(-1, 1, .05)){
  lines(seq(0, 100, .05), prev_sv_simetric(seq(0, 100, .05), 100, lambdas), 
        col=ifelse(lambdas>0, 1, 2))
}

# get lambda to get these prev
get_best_lambda = function(lambda, S, alpha){
  Area_lambda = integrate(f = prev_sv_simetric, lower =  0, upper = S, subdivisions = 20000, 
                          y = S, lambda = lambda)$value
  (Area_lambda-alpha*S)^2
}

# get prevalence of this guy
prev_lambda_alpha_y <- function( x = seq(0,omega,by=delta), y = max(x), 
                                  alpha = .4, S = 20, omega = 110,delta = .05){
  # target prev
  lambda = optimize(get_best_lambda, lower=-2, upper=2, alpha = alpha, S = S, tol = 1e-12)$minimum 
  prev_S = prev_sv_simetric(seq(0,S,delta), y = S, lambda = lambda)
  prev_x = data.frame(x = x, prev = 0)
  prev_x$prev[prev_x$x>=(y-S) & prev_x$x<(y+delta)] =  prev_S
  prev_x$prev[prev_x$x > y] <- NA
  return(list(prev_x = prev_x, lambda = lambda))
  }

plot(NULL, type = "n", xlim = c(0,50), ylim = c(0,1))
lines(prev_lambda_alpha_y(y=50,alpha=.8,S=20)$prev_x)

# get prevalences for all
prev_x_from_lambda_alpha <- function( lives_lx, 
                                      x = seq(0, omega, by = delta), 
                                      alphas = .4, 
                                      Ss = 20, 
                                      delta = .05, 
                                      omega = 110){
  n <- length(lives_lx)
  if (length(alphas) == 1){
    alphas <- rep(alphas, n)
  }
  if (length(Ss) == 1){
    Ss <- rep(Ss, n)
  }
  # get prevalence trajectory for each lifespan
  Prev_lifespan <- matrix(0, ncol = length(x), nrow = length(lives_lx))
  for (i in 1:length(lives_lx)){
    # i = 1
    Prev_lifespan[i,] <- prev_lambda_alpha_y( x = x, 
                                              y = lives_lx[i], 
                                              alpha = alphas[i], 
                                              S = Ss[i],
                                              delta = delta,
                                              omega = omega)$prev_x$prev
  }
  Prev_x <- colMeans(Prev_lifespan, na.rm = TRUE)
  return(list(Prev_x = Prev_x, Prev_lifespan = Prev_lifespan, x = x))
}

# get 
bin_prev_x <- function(prev_x, 
                       xin = seq(0,110,by=.05), 
                       xout = 0:100){
  # prev_x = Prev_lifespan
  delta_implied = diff(xin)[1]
  xout_index = trunc(xin)
  prev_xout = matrix(0, nrow=nrow(prev_x), ncol=length(xout))
  for (ls in 1:nrow(prev_x)){
    #ls=1
    age_d = sum(!is.na(prev_x[ls,])) * delta_implied 
    for (age in xout){
      # age = 16
      # have to take account for deaths between integer ages
      prev_xout[ls, age+1] = sum(prev_x[ls, xout_index==age], na.rm = T) * delta_implied #/ min(1,age_d-age, na.rm = T)
    }
  }
  prev_xout
}

# sum vertical
prev_cross_x = function(Prev_lifespan_bin, lives_lx, xout){
  alives = sapply(xout, function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
  Prev_x = pmin(colSums(Prev_lifespan_bin, na.rm = T)/alives, 1)
  Prev_x
}

# together
prev_x_final = function(lives_lx, x, alphas, Ss, delta, omega, xout){
  # get prev matrix
  Prev_lifespan = prev_x_from_lambda_alpha(lives_lx, x, alphas, Ss, delta, omega)$Prev_lifespan
  # group by integer ages
  Prev_lifespan_bin = bin_prev_x(Prev_lifespan, xin = x, xout = xout)
  # person-year survivors in each integer age (decimal y gives the portion also)
  Prev_x = prev_cross_x(Prev_lifespan_bin, lives_lx, xout)
  Prev_x
}

###---------------------- optimization given a prevalence target

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
omega <- 110; delta <- .05
lives_lx <- life_bins(lx,x=0:110,q)
lives_lx <- round(lives_lx * 1/delta) * delta
lives_lx[lives_lx > omega] <- omega
S = 20 #same S for every
x = seq(0, omega, by = delta)
xout = 0:100
Pi_Obs = .8/(1+exp(-.1*(xout-70))) # works for lambdas > .05 !
# plot(xout, Pi_Obs, ylim=c(0,1))


### ojete function
fun_opt_both = function(par, 
                        x , lx, S,
                        delta,  omega , q , xout, xobj,
                        Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:omega, probs = q)[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  # alphas = pmin(my_poly(pars[1:2], lives_lx), 2)
  alphas = pmin(my_poly(c(0,par), lives_lx), 2)
  # vectorize in case
  Ss = rep(S, length.out = length(lives_lx))
  
  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

optim.s = optimize(    # starting pars
                    f = fun_opt_both,     # minimize this
                    lower = 0.0001,
                    upper = .01,
                    S = S, # doy una constante
                    lx = lx, # mortality 
                    q = seq(1,0,by=-.01), #seq(1,0,by=-.01), 
                    x = seq(0, omega, by = delta), 
                    delta = .05,  omega = 110, 
                    xout = xout, # integer ages (could be implicit in future)
                    xobj = 30:90, # fit priority
                    Pi_obs = Pi_Obs)

alphas_optim = pmin(my_poly(c(0,optim.s$minimum), lives_lx), 2) #; plot(lives_lx, alphas_optim)
Prev_hat = prev_x_final(lives_lx, x, 
                        Ss = rep(S, length.out = length(lives_lx)), 
                        alphas = alphas_optim,
                        delta, omega, xout)
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
points(xout, Prev_hat, t='l', col = 2)

