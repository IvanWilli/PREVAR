### IW: Baaad performance

# set age breaks
Pi_Obs = .8/(1+exp(-.1*(xout-70)))
plot(Pi_Obs)
plot(diff(diff(Pi_Obs)))
break_age = which(diff(diff(Pi_Obs))==max(diff(diff(Pi_Obs))))+2
break_age = round(break_age/5)*5
break_ages = c(0, 30, 60, 85, 100)

# obj func
fun_opt_both = function(pars, 
                        x , lx, S,
                        delta,  omega , q , xout, 
                        xobj,
                        Pi_obs){
  # lives
  lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
  lives_lx <- round(lives_lx * (1/delta))*delta
  lives_lx[lives_lx>omega] <- omega
  # pars = param_0
  alphas = pmin(my_poly(pars[1:2], lives_lx[lives_lx >= min(xobj) & lives_lx <= max(xobj)]), 2)
  lambdas = pmin(my_poly(pars[3:4], lives_lx[lives_lx >= min(xobj) & lives_lx <= max(xobj)]), 7)
  Ss = rep(S, length.out = length(lives_lx[lives_lx >= min(xobj) & lives_lx <= max(xobj)]))
  Prev_hat = prev_x_final(lives_lx[lives_lx >= min(xobj) & lives_lx <= max(xobj)], x, lambdas, alphas, Ss, delta, omega, xout)
  sum((Pi_obs[xobj+1] - Prev_hat[xobj+1]) ^2, na.rm = T)
}


param_Lbounds = c(0.0001,  .00001,    .000001,   .0000001)
param_Ubounds = c(1,        .1,       5,         .05)
S = 20
xout = 0:100

Optims_ranges = data.frame(agei=NA, agef=NA, alpha1=NA, alpha2=NA, lambda1=NA, lambda2=NA, objfun = NA)
set.seed(45)
for (i in 1:(length(break_ages)-1)){
  # i = 3
  param_0 =  param_Lbounds + (param_Ubounds-param_Lbounds) * rep(runif(1),4)
  xobj = break_ages[i]:break_ages[i+1]
  optim_i = solnp(pars = param_0,    # starting pars
                  fun = fun_opt_both,     # minimize this
                  LB = param_Lbounds, 
                  UB = param_Ubounds,
                  S = S,
                  lx = lx, 
                  q = seq(1,0,by=-.01), 
                  x = seq(0, omega, by = delta), 
                  delta = .05,  omega = 110, 
                  xout = xout, 
                  xobj = xobj, 
                  Pi_obs = Pi_Obs,
                  control = list(delta = .01, rho = 1.5, outer.iter = 50))
  Optims_ranges = rbind(Optims_ranges,
                         data.frame(agei=break_ages[i], agef=break_ages[i+1], 
                                    alpha1=optim_i$pars[1], alpha2=optim_i$pars[2], 
                                    lambda1=optim_i$pars[3], lambda2=optim_i$pars[4], 
                                    objfun = optim_i$values[1]))
}


# Fit

Optims_ranges = Optims_ranges[!is.na(Optims_ranges$agei),]
Prev_hat = matrix(0, ncol=length(xout), nrow=nrow(Optims_ranges))
for (j in 1:nrow(Optims_ranges)){
  # j = 2
  alphas_optim = pmin(my_poly(Optims_ranges[j, 3:4], lives_lx), 2)
  lambdas_optim = pmin(my_poly(Optims_ranges[j, 5:6], lives_lx), 2)
  Ss = rep(S, length.out = length(lives_lx))
  Prev_hat[j,] = prev_x_final(lives_lx, x, lambdas_optim, alphas_optim, Ss, delta, omega, xout)
}

plot(0:100, Pi_Obs, xlim = c(0,100), ylim = c(0,1), xlab = 'Ages', ylab = 'Prevalence', main = 'Scenarios Fits')
points(xout, Prev_hat[1,], t='l', col = 2)
points(xout, Prev_hat[2,], t='l', col = 3)
points(xout, Prev_hat[3,], t='l', col = 4)
points(xout, Prev_hat[4,], t='l', col = 5)

### IW: I think the starting parameter sensitivity makes that older segments needs greater starting points
### maybe too much local minimums in the combination lambda-alpha surface