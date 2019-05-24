############ what is the 'disability mass' by time_to_death to replicate prev #########################
############ get lambda and alpha ###############

# objective function
fun_opt_both = function(pars, 
                          family = 'lineal', 
                          x , lx, S,
                          delta,  omega , q , xout, xobj,
                          Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  a_alpha = pars[1]
  b_alpha = pars[2]
  a_lambda = pars[3]
  b_lambda = pars[4]
  alphas = pmin(alpha_fun(x = lives_lx, a = a_alpha, b = b_alpha, family = family), 2)
  lambdas = pmin(lambda_fun(x = lives_lx, a = a_lambda, b = b_lambda, family = family), 10)
  
  # vectorize in case
  Ss = rep(S, length.out = length(lives_lx))
  
  # get Prev_hat
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  
  # obj funct en rango objetivo
  sum((Pi_obs[xobj+1] -Prev_hat[xobj+1]) ^2)
}

# Set family and bounds

family='lineal'
if (family=='lineal'){
  param_0 =       c(.5,      .01,    0.00001,    .01)
  param_Lbounds = c(0.0001,  .00001, .000001,    .000001)
  param_Ubounds = c(2,       .05,    1,        .05  )
}
S = 20 #same S for every
xout = 0:100
xobj = 30:80 # fit priority
# Possible observed prevalences
      Pi_Obs = .6/(1+exp(-.1*(xout-80)))
      # Pi_Obs = .005*exp(.05*xout)
      # plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1))
        
# Optimization
# require(Rsolnp)
optim.s = solnp(param_0,    # starting pars
                fun = fun_opt_both,     # minimize this
                LB = param_Lbounds, 
                UB = param_Ubounds,
                family = family,
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
alphas_optim = alpha_fun(x = lives_lx, a = optim.s$pars[1], b = optim.s$pars[2], family = family)
lambdas_optim = lambda_fun(x = lives_lx, a = optim.s$pars[3], b = optim.s$pars[4], family = family)
Prev_hat = prev_x_final(lives_lx, x, 
                        lambdas = lambdas_optim, 
                        Ss = rep(S, length.out = length(lives_lx)), 
                        alphas = alphas_optim,
                        delta, omega, xout)
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
points(xout, Prev_hat, t='l', col = 2)

# IW = Good fit. Only tried with linealin both lambda and alpha