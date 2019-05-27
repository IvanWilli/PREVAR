############ what is the 'disability mass' by time_to_death to replicate prev #########################

############ get lambda and alpha ###############
me <- system("hostname",intern=TRUE) 
if (me == "tim-ThinkPad-L440"){
  setwd("/home/tim/git/PREVAR/PREVAR")
}
#install.packages("Rsolnp")
library(Rsolnp)

source(file.path("R","TTDFunctions.R"))
# objective function
# objective function
fun_opt_both = function(pars, 
                        n_poly = c(1,1), 
                        x , lx, S,
                        delta,  omega , q , xout, xobj,
                        Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  alphas = pmin(my_poly(pars[1:(n_poly[1]+1)], lives_lx), 2)
  lambdas = pmin(my_poly(pars[3:(n_poly[2]+1)], lives_lx), 2)
  
  # vectorize in case
  Ss = rep(S, length.out = length(lives_lx))
  
  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

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
# Set family and bounds

q        <- seq(1,0,by=-.005)
lives_lx <- life_bins(lx,x=0:110,q)
lives_lx <- round(lives_lx * 1/delta) * delta
lives_lx[lives_lx > omega] <- omega

# just for lineal case
param_0 =       c(.5,      .01,    0.00001,    .01)
param_Lbounds = c(0.0001,  .00001, .000001,    .000001)
param_Ubounds = c(2,       .05,    1,        .05  )

S = 20 #same S for every
xout = 0:100
xobj = 30:80 # fit priority
# Possible observed prevalences
Pi_Obs = .8/(1+exp(-.1*(xout-80)))
# Pi_Obs = .005*exp(.05*xout)
# plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1))

# Optimization
# require(Rsolnp)
optim.s = solnp(param_0,    # starting pars
                fun = fun_opt_both,     # minimize this
                LB = param_Lbounds, 
                UB = param_Ubounds,
                n_poly = c(1,1),
                S = S, # doy una constante
                lx = lx, # mortality 
                q = seq(1,0,by=-.01), 
                x = seq(0, omega, by = delta), 
                delta = .05,  omega = 110, 
                xout = xout, # integer ages (could be implicit in future)
                xobj = xobj, # fit priority
                Pi_obs = Pi_Obs,
                control = list(delta = .01, rho = 1.5, outer.iter = 50))

# Visual Goodness of Fit
alphas_optim = pmin(my_poly(optim.s$pars[1:(n_poly[1]+1)], lives_lx), 2)
lambdas_optim = pmin(my_poly(optim.s$pars[3:(n_poly[2]+1)], lives_lx), 2)
Prev_hat = prev_x_final(lives_lx, x, 
                        lambdas = lambdas_optim, 
                        Ss = rep(S, length.out = length(lives_lx)), 
                        alphas = alphas_optim,
                        delta, omega, xout)
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
points(xout, Prev_hat, t='l', col = 2)

# IW = Good fit. Only tried with linealin both lambda and alpha
# ONLY LINEAL, 















fun_opt_poly = function(pars, 
                        x, 
                        lx, 
                        lambda = .5,
                        S = 30,
                        delta = .05, 
                        omega = 110, 
                        q = seq(1,0,by=.005) , 
                        xout = 0:110, 
                        xobj = 30:80,
                        Pi_obs){
  # lives
  lives_lx <- life_bins(lx, x, probs = q )[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  
  
  # vectorize in case
  if (length(S) == 1){
    Ss = rep(S, length.out = length(lives_lx))
  }
  
  alphas   <- alpha_poly(lives_lx,pars[1:3])
  #lambdas  <- lambda_poly(lives_lx, pars[4:6])
  lambdas <- rep(lambda,length(lives_lx))
  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

opt <- optim(c(0,0,0),
             fun_opt_poly,
             lower=-.5,
             upper=.5,
             x=0:110,
             lx=lx,
             S=30,
             lambda=.2,
             delta=.05,
             omega=110,
             q = seq(1,0,by=-.005),
             xout = 0:110,
             xobj = 30:80,
             Pi_obs = Pi_Obs)
str(opt)
plot(alpha_poly(0:110,opt$par[1:3]))
#plot(lambda_poly(0:110,opt$par[4:6]))


Prev_hat = prev_x_final(lives_lx, x, 
                        lambdas = rep(.5,length(lives_lx)), 
                        Ss = rep(30, length.out = length(lives_lx)), 
                        alphas = alpha_poly(lives_lx,opt$par[1:3]),
                        delta, omega, xout)
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
points(xout, Prev_hat, t='l', col = 2)