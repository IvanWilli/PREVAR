fun_opt = function(pars, x , lx, 
                   delta,  omega , q , xout, xobj,
                   Pi_obs){

          # lives
          lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
          lives_lx[lives_lx>omega] <- omega
          
          # pars
          # pars = c(1, .05, .1, 0.01); S = 15
          # pars = optim.s$pars
          al = pars[1]
          bl = pars[2]
          aa = pars[3]
          ba = pars[4] 
          S = pars[5] 
          lambdas = lives_lx * bl + al
          alphas = lives_lx * ba + aa
          Ss = rep(S, length.out = length(lives_lx))
          # plot(lives_lx, lambdas, xlim=c(0,100), ylim=c(0,10)); plot(lives_lx, alphas, xlim=c(0,100), ylim=c(0,2))
          
          # get prev matrix
          Prev_lifespan = prev_x_from_lambda_alpha(lives_lx, x, lambdas, alphas, Ss, delta, omega)$Prev_lifespan

          # group by integer ages
          Prev_lifespan_bin = bin_prev_x(Prev_lifespan, xin = x, xout = xout)
          
          # person-year survivors in each integer age (decimal y gives the portion also)
          Prev_x = prev_cross_x(Prev_lifespan_bin, lives_lx, xout)
          # plot(xobj, Pi_obs[xobj+1], xlim=range(xobj))
          # lines(xobj, Prev_x[xobj+1], xlim=range(xobj))
          
          # obj funct en rango objetivo
          sum((Pi_obs[xobj+1] -Prev_x[xobj+1]) ^2)
          # if(opti_result=='opti') {
          #   result = sum((Pi_obs[xobj+1] -Prev_x[xobj+1]) ^2)}
          # if(opti_result=='prev') {
          #   result = Prev_x[xobj+1]}
          # return(result)
}


param_Lbounds = c(0,   0.001,   0.0001,  .001,  5)
param_Ubounds = c(5,   .1,      2,       .1,   20)
param_0 = c(1, .05, .1, 0.01, 10)

# require(Rsolnp)
optim.s = solnp(param_0,    # starting pars
                fun = fun_opt,     # minimize this
                LB = param_Lbounds, 
                UB = param_Ubounds,
                x = seq(0, omega, by = delta), lx = lx, 
                delta = .05,  omega = 110, q = seq(1,0,by=-.01), 
                xout=0:100, xobj = 30:80,
                Pi_obs =.4/(1+exp(-.1*(xout-80))),
                control = list(delta = .1, rho = 1.5, outer.iter = 50)
                )









