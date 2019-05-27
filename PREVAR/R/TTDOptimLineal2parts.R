# set age breaks
plot(diff(diff(Pi_Obs)))
break_age = which(diff(diff(Pi_Obs))==max(diff(diff(Pi_Obs))))+2
break_age = round(break_age/5)*5

# obj func
fun_opt_both = function(pars, 
                        x , lx, S,
                        delta,  omega , q , xout, xobj,
                        Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  
  # pars = param_0
  alphas = pmin(my_poly(pars[1:2], lives_lx), 2)
  lambdas = pmin(my_poly(pars[3:4], lives_lx), 7)
  Ss = rep(S, length.out = length(lives_lx))
  Prev_hat = prev_x_final(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout)
  sum((Pi_obs[range(xobj)+1] - Prev_hat[range(xobj)+1]) ^2)
}

# set bounds
param_0 =      c(.5,        .01,      1,         .00000001)
param_Lbounds = c(0.0001,  .00001,    .000001,   -.05)
param_Ubounds = c(2,        .1,       5,         .05)
Pi_Obs = .6/(1+exp(-.1*(xout-80)))
S = 20
xout = 0:100
xobj = 30:90

# 2 optims
optim.s = solnp(pars = param_0,    # starting pars
                fun = fun_opt_both,     # minimize this
                LB = param_Lbounds, 
                UB = param_Ubounds,
                S = S, # doy una constante
                lx = lx, # mortality 
                q = seq(1,0,by=-.01), 
                x = seq(0, omega, by = delta), 
                delta = .05,  omega = 110, 
                xout = xout, # integer ages (could be implicit in future)
                xobj = min(xobj):break_age, # fit priority
                Pi_obs = Pi_Obs,
                control = list(delta = .01, rho = 1.5, outer.iter = 50))
optim.s2 = solnp(pars = param_0,    # starting pars
                fun = fun_opt_both,     # minimize this
                LB = param_Lbounds, 
                UB = param_Ubounds,
                S = S, # doy una constante
                lx = lx, # mortality 
                q = seq(1,0,by=-.01), 
                x = seq(0, omega, by = delta), 
                delta = .05,  omega = 110, 
                xout = xout, # integer ages (could be implicit in future)
                xobj = (break_age+1):max(xobj), # fit priority
                Pi_obs = Pi_Obs,
                control = list(delta = .01, rho = 1.5, outer.iter = 50))

# Fit
alphas1_optim = pmin(my_poly(optim.s$pars[1:2], lives_lx), 2)
lambdas1_optim = pmin(my_poly(optim.s$pars[3:4], lives_lx), 2)
alphas2_optim = pmin(my_poly(optim.s2$pars[1:2], lives_lx), 2)
lambdas2_optim = pmin(my_poly(optim.s2$pars[3:4], lives_lx), 2)
Ss = rep(S, length.out = length(lives_lx))
Prev_hat1 = prev_x_final(lives_lx, x, lambdas1_optim, alphas1_optim, Ss, delta, omega, xout)
Prev_hat2 = prev_x_final(lives_lx, x, lambdas2_optim, alphas2_optim, Ss, delta, omega, xout)
Prev_hat = c(Prev_hat1[min(xobj):break_age+1], Prev_hat2[(break_age+1):max(xobj)+1])
plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
abline(v = range(xobj), col = 'grey', lty = 2)
points(xobj, Prev_hat, t='l', col = 2)
