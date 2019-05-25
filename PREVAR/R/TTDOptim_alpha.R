############ Given levels of lambda, what is the 'disability mass' by time_to_death to replicate prev #########################

# objective function
fun_opt_lambda = function(pars, 
                          lambdas, 
                          family = 'lineal', x , lx, S,
                          delta,  omega , q , xout, xobj,
                          Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  a_alpha = pars[1]
  b_alpha = pars[2]
  alphas = pmin(alpha_fun(x = lives_lx, a = a_alpha, b = b_alpha, family = family), 2)
  
  # vectorize in case
  lambdas = rep(lambdas, length.out = length(lives_lx))
  Ss = rep(S, length.out = length(lives_lx))

  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

# Diferent scenarios of lambda and family

lambdas_scenarios = list(0, 
                         1, 
                         seq(5, 0, -.05)[1:length(lives_lx)])
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
    lambda = lambdas_scenarios[[i]] # could be function by age
    xout = 0:100
    Pi_Obs = .6/(1+exp(-.1*(xout-80))) # observed prev
    
    # Optimization
    # require(Rsolnp)
    optim.s = solnp(param_0,    # starting pars
                    fun = fun_opt_lambda,     # minimize this
                    LB = param_Lbounds, 
                    UB = param_Ubounds,
                    family = family,
                    lambdas = lambda, S = S, # doy una constante
                    lx = lx, # mortality 
                    q = seq(1,0,by=-.01), 
                    x = seq(0, omega, by = delta), 
                    delta = .05,  omega = 110, 
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
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
cols = c(4,2,4,2,4,2); ltys = c(1,1,2,2,3,3)
for (i in 1:nrow(Fits_lambdas_scenarios)){
  # i=6
  alphas_optim = alpha_fun(x = lives_lx, a = Fits_lambdas_scenarios$a_alpha[i], 
                           b = Fits_lambdas_scenarios$b_alpha[i], 
                           family = Fits_lambdas_scenarios$family[i])
  Prev_hat = prev_x_final(lives_lx, x, 
                          lambdas = Fits_lambdas_scenarios$lambdas_esc[i], 
                          Ss = rep(S, length.out = length(lives_lx)), 
                          alphas = alphas_optim,
                          delta, omega, xout)
  points(xout, Prev_hat, t='l', col = cols[i], lty = ltys[i])
}
legend('topleft', paste(Fits_lambdas_scenarios$family,Fits_lambdas_scenarios$lambdas_esc, sep=' - esc Lambda '),  
       col=c(4,2,4,2,4,2), bty='n', lty = c(1,1,2,2,3,3))
  
# IW: Lineal fit seems good for all lambda scenarios
#     Exponential family doesn't fit because of the shape. We'll have to try quadratic or another
#     Prevalence slowing down 80+ dificult to catch
#     iteration is too slow, maybe 'integrate' function must be replaced



