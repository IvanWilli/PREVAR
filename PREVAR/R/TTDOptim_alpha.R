############ Given levels of lambda, what is the 'disability mass' by time_to_death to replicate prev #########################

me <- system("hostname",intern=TRUE) 
if (me == "tim-ThinkPad-L440"){
	setwd("/home/tim/git/PREVAR/PREVAR")
}
#install.packages("Rsolnp")
library(Rsolnp)

source(file.path("R","TTDFunctions.R"))
# objective function
fun_opt_lambda = function(pars, 
                          lambdas, 
                          family = 'lineal', x , lx, S,
                          delta,  omega , q , xout, xobj,
                          Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  a_alpha = pars[1]
  b_alpha = pars[2]
  alphas = pmin(alpha_fun(x = lives_lx, a = a_alpha, b = b_alpha, family = family), 2)
  
  # vectorize in case
  lambdas = rep(lambdas, length.out = length(lives_lx))
  if (length(S) == 1){
	  Ss = rep(S, length.out = length(lives_lx))
  } else {
	  Ss <- S
  }

  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

# Diferent scenarios of lambda and family

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

q        <- seq(1,0,by=-.005)
lives_lx <- life_bins(lx,x=0:110,q)

delta    <- 0.05
lives_lx <- round(lives_lx * 1/delta) * delta
omega    <-  110
lives_lx[lives_lx > omega] <- omega
x <- seq(0,omega,by=delta)
# no e0 penalty for small delta
#mean(life_bins(lx,x=0:110,q))-mean(lives_lx)

lambdas_scenarios = list(0, 
                         1, 
                         seq(5, .5, length=length(lives_lx)))
families = c('lineal', 'exp')

# Results matrix
Fits_lambdas_scenarios = data.frame(family = 0, lambdas_esc = 0, a_alpha = 0, b_alpha = 0)

for (i in 1:length(lambdas_scenarios)) {
  
  for (family in families){
    
    # set family, starting pars and bounds (depends of family)
    if (family=='lineal'){
      param_0 =       c(.5,      .01)
      param_Lbounds = c(0.0001,  .00001)
      param_Ubounds = c(2,       .05 )
    }
    if (family=='exp'){
      param_0 =       c(.01,      .05)
      param_Lbounds = c(0.0001,  .00001)
      param_Ubounds = c(2,       .07 )
    }
    # plot(0:100, alpha_fun(x = 0:100, a = param_0[1], b = param_0[2], family = family), xlim=c(0,100), ylim=c(0,2))
    S = 20 #same S for every
	# Wondering now if things would work better with a nice big 40 instead? This would allow small lambdas
    # to have more effect? hmm.
    lambda = lambdas_scenarios[[i]] # could be function by age
    xout = 0:110
    Pi_Obs = .6/(1+exp(-.1*(xout-80))) # observed prev
    
    # Optimization
    # require(Rsolnp)
    optim.s = solnp(param_0,    # starting pars
                    fun = fun_opt_lambda,     # minimize this
                    LB = param_Lbounds, 
                    UB = param_Ubounds,
                    family = family,
                    lambdas = lambda, 
					S = 40, # doy una constante
                    lx = lx, # mortality 
                    q = seq(1,0,by=-.01), 
                    delta = .05,  
					omega = 110, 
					x = seq(0, omega, by = delta), 
                    xout = xout, # integer ages (could be implicit in future)
                    xobj = 30:80, # fit priority
                    Pi_obs = Pi_Obs,
                    control = list(delta = .001, rho = 1.5, outer.iter = 50)
    )
    Fits_lambdas_scenarios = rbind(Fits_lambdas_scenarios, 
                                   data.frame(family = family, lambdas_esc = i, 
                                              a_alpha = optim.s$pars[1], b_alpha = optim.s$pars[2]))
  }
  if (i == length(lambdas_scenarios)) {
    Fits_lambdas_scenarios = Fits_lambdas_scenarios[-1,]
  }
}
# takes 7 minutes

# Visual Goodness of Fit
plot(xout, Pi_Obs, xlim = c(0,100), ylim = c(0,1), 
		xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits', 
		type = 'l', lwd = 2)
abline(v = c(30,80), col = 'grey', lty = 2)
cols = c(4,2,4,2,4,2); ltys = c(1,1,2,2,3,3)
for (i in 1:nrow(Fits_lambdas_scenarios)){
  # i=6
  alphas_optim = alpha_fun(x = lives_lx, a = Fits_lambdas_scenarios$a_alpha[i], 
                           b = Fits_lambdas_scenarios$b_alpha[i], 
                           family = Fits_lambdas_scenarios$family[i])
  Prev_hat = prev_x_final(lives_lx, x, 
                          lambdas = Fits_lambdas_scenarios$lambdas_esc[i], 
                          Ss = rep(40, length.out = length(lives_lx)), 
                          alphas = alphas_optim,
                          delta=0.05, omega=110, xout=0:110)
  lines(0:110, Prev_hat, col = cols[i], lty = ltys[i])
}
legend('topleft', paste(Fits_lambdas_scenarios$family,Fits_lambdas_scenarios$lambdas_esc, sep=' - esc Lambda '),  
       col=c(4,2,4,2,4,2), bty='n', lty = c(1,1,2,2,3,3))
  
# IW: Lineal fit seems good for all lambda scenarios
#     Exponential family doesn't fit because of the shape. We'll have to try quadratic or another
#     Prevalence slowing down 80+ dificult to catch
#     iteration is too slow, maybe 'integrate' function must be replaced
# TR: replaced integrate() with analytic prev_sv_int() function 
#     



