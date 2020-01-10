# health & mortality INPUTS
lx_orig <- c(100000L, 99326L, 99281L, 99249L, 99226L, 99208L, 99194L, 99182L, 
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

Pi_Obs1 = 1/(1+exp(-.05*(0:110-140)))
Pi_Obs4 = 1/(1+.002*exp(-.206*(0:110-140))) + .01
plot(Pi_Obs1, t="l");lines(Pi_Obs4,lty=2); abline(v=c(50,100),lty=3)
Prev = data.frame(x=0:110, Prev = Pi_Obs4)
lt = data.frame(x=0:110, l=lx_orig)
lt = lt[lt$x>=50,] 
lx = lt$l; x = lt$x
Prev = Prev[Prev$x>=50,"Prev"] 

# some runs
ttd_run <- ttd_ineq(x = x, lx = lx, bin = 1, Prev, S = 5, 
               lambda_limit = c(0,100), ages_care = 50:99)
ttd_run <- ttd_ineq(x = x, lx = lx, bin = 1, Prev, S = 10, 
                    lambda_limit = c(0,100), ages_care = 50:99)
ttd_run <- ttd_ineq(x = x, lx = lx, bin = .5, Prev, S = 5, 
                    lambda_limit = c(0,100), ages_care = 50:99)


# function
ttd_ineq <- function(x, lx, bin = 1, Prev, S = 10, 
                     lambda_limit = c(0,100), ages_care = 50:99){
  
  # set
  x_int <- x
  x     <- seq(x[1],x[length(x)],bin)
  lx    <- splinefun(x_int, lx, method = "monoH.FC")(x)
  n     <- 1:(length(lx))
  dx    <- c(-diff(lx),lx[length(lx)])
  Prev  <- splinefun(x = x_int, y = Prev)(x)
  x     <- c(x,  x[length(x)]+bin)
  lx    <- c(lx,0)
  
  # lambdas & matrix prev
  matrix_prev_YP <- matrix(0, length(x)-1, length(lx)-1)
  lambdas <- c()
  for(y in rev(n)){
    # y = 600
    lambdas[y]  <- optimize(get_lambda, 
                            S = S, 
                            lower = lambda_limit[1], upper = lambda_limit[2], 
                            x_ini = x[y], x_fin = x[y+1],
                            PrevObs_x = Prev[y] * lx[y] - sum(matrix_prev_YP[,y]), 
                            lives = dx[y])$minimum
    
    matrix_prev_YP[y, y:max(y+1-S/bin,1)] <- sapply(seq(x[y],max(x[y]-S+bin,x[1]), -bin), 
                                                     FUN = function(s) 
                                                     integrate(f = prev_sv_simetric, 
                                                       lower = s, upper = (s+bin),  
                                                       y = x[y+1], S = S, 
                                                       lambda = lambdas[y],
                                                       subdivisions = 20000)$value  * dx[y]
                                                    ) 
  }
  
  # first fit
  Prev_fit  <- colSums(matrix_prev_YP[,-ncol(matrix_prev_YP)])/lx[-length(lx)]
  Prev_fit <- colSums(matrix_prev_YP)/lx[-length(lx)]
  
  
  # beta fit of lambdas
  pars_fit <- optim(par = pars_init, beta_min_lala, ages = x_int, ages_care = ages_care, 
                    lambda_obs = lambdas)$par
  x_dist   <- seq(0, 1, length = length(lambdas))
  lambdas_fit <- pars_fit["scale"] * dbeta(x_dist, shape1 =  pars_fit["shape1"], 
                                                   shape2 =  pars_fit["shape2"])
  
  
  matrix_prev_YP_hat <- matrix(0, length(x)-1, length(lx)-1)
  for(y in rev(n)){
    matrix_prev_YP_hat[y, y:max(y+1-S/bin,1)] <- sapply(seq(x[y],max(x[y]-S+bin,x[1]), -bin), 
                                                    FUN = function(s) 
                                                    integrate(f = prev_sv_simetric, 
                                                        lower = s, upper = (s+bin),  
                                                        y = x[y+1], S = S, 
                                                        lambda = lambdas_fit[y],
                                                        subdivisions = 20000)$value * dx[y])  
  }
  
  # graphs
  x = x[-length(x)] 
  Prev_hat  <- colSums(matrix_prev_YP_hat) / lx[-length(lx)]
  plot(x, Prev, pch=15, cex=.5);lines(x, Prev_fit);lines(x, Prev_hat,col=2); abline(v = range(ages_care), lty=2) 
  legend("topleft", c("Prev Obs","Prev fit lambdas", "Prev fit betas", "Ages that care"),
        lty=c(NA,1,1,2), col=c(1,1,2,1), pch=c(15,NA,NA,NA), bty="n")
  prev_graph <- recordPlot()
  plot(x, lambdas); lines(x, lambdas_fit)
  legend("topleft", c("lambdas", "lambdas fitted by beta dist"), lty=c(1,NA), pch=c(NA,1), cex=.8, bty="n")
  lambda_graph <- recordPlot()
  graphs = list(prev_fit = prev_graph, lambda_fit = lambda_graph)
  
  # goodness of fit in ages care
  GoF       <- list(MSE_lambda = sum((Prev[ages_care %in% x]-Prev_fit[ages_care %in% x])^2/Prev[ages_care %in% x]),
                    MSE_betas  = sum((Prev[ages_care %in% x]-Prev_hat[ages_care %in% x])^2/Prev[ages_care %in% x]))

  # ineq measures
  variance = var(rowSums(matrix_prev_YP))
  entropy_lt = sum(-log(dx/lx[1])*dx/lx[1]) 
  entropy_lt_prev = sum(rowSums(-log(matrix_prev_YP) * matrix_prev_YP, na.rm=T) * dx/lx[1])
  ineq_index = list(variance = variance, 
                    entropy_lt=entropy_lt, entropy_lt_prev=entropy_lt_prev, 
                    entropy = entropy_lt+entropy_lt_prev)
  
  return(list(matrix_prev_YP = matrix_prev_YP, 
              Prev_fit = Prev_fit, Prev_hat = Prev_hat, 
              GoF = GoF, 
              ineq_index = ineq_index, 
              graphs = graphs))
}

# obj function
beta_min_lala <- function(pars = c(scale = 1, shape1 = 10, shape2 = 5), 
                          lambda_obs, ages, ages_care){
  
  x_dist   <- seq(0, 1, length = length(lambda_obs))
  x_care_dist <- seq((min(ages_care)-min(ages))/diff(range(ages)),
                     (max(ages_care)-min(ages))/diff(range(ages)),
                     by = diff(x_dist)[1])
  lambdas <- pars["scale"] * 
    dbeta(x_care_dist, 
          shape1 =  pars["shape1"], shape2 =  pars["shape2"])
  sum((lambda_obs[x_dist %in% x_care_dist] - lambdas[x_dist %in% x_care_dist])^2)
}

# beta fit
pars_fit <- sapply(1:ncol(lambdas_adj), 
                   function(s) optim(pars_init, 
                                     beta_min_lala,
                                     ages, ages_care,
                                     lambda_obs = lambdas_adj[,s]$lambda)$par)
