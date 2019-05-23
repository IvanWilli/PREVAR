# Prev function
# x decimal evaluated age
# y decimal death age
# lambda is compression parameter | higher = more compressed
prev_sv = function(x, y, lambda=.3){ # x vector of ages at death
  p = ifelse(is.infinite(x*exp(x*lambda)), 0, x*exp(x*lambda))/(y*exp(y*lambda))
  p = ifelse(p>1,0,p)
  p
}

# interpolate age at survival
life_bins <- function(lx, x, probs =  seq(0,1,by=.05)){
    lx <- lx / lx[1]
    splinefun(x~lx, method = "monoH.FC")(probs)
}
    
    
# a utility to help optimize for t_star, which scales the prevalence area by 
# shifting the prevalence curve (implied by alpha and lambda) left or right out of the window S
get_shift = function(t_star, alpha, lambda, S){ 
      
      # https://www.wolframalpha.com/input/?i=%5Cint_0%5Et%7Bx%2Ft*%5Cexp((x-t)*%5Clambda)%7D
      # Area_hat = (lambda * t_star + exp(lambda * -t_star) - 1) / (lambda^2) / t_star
      
      if (sign(t_star) == 1){
        Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S, lambda = lambda)$value + t_star
        
      } else {
        if (sign(t_star) == -1){
          Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S + t_star, lambda = lambda)$value
        } 
      }
      
      ((alpha * (S/2)) - (Area_lambda_S)) ^ 2
    }


prev_lambda_alpha_y <- function(
                                x = seq(0,omega,by=delta), 
                                y = max(x), 
                                lambda = 0, 
                                alpha = 1.5, 
                                S = 20, 
                                omega = 110,
                                delta = .05){
  
  # start with....
  # 0      y-S      y    t*
  # |______|________|______|
  #            S      t*-S
  
  t_star = optimize(get_shift, interval = c(-S, S), alpha = alpha, lambda = lambda, S = S, tol = 1e-12)$minimum
  
  t_star <- round(t_star * 1/delta) * delta
  
  prev   <- x * 0
  if (sign(t_star) == 1){
    prev_x <- prev_sv(x = seq(t_star, S, by = delta), y = S, lambda = lambda) 
    x_implied <- seq((y - S), (y - t_star), by = delta)
    prev[x >= (y - S) & x <= (y - t_star)] <- prev_x[x_implied>=0]
    prev[x > (y - t_star) & x <= y] <- 1
  }
  if (sign(t_star) == -1){
    prev_x <- prev_sv(x = seq(0, S + t_star, by = delta), y = S, lambda = lambda)
    prev[x >= (y-S-t_star) & x <= (y)] <- prev_x
  }
  if (sign(t_star) == 0){
    prev_x <- prev_sv(x = seq(0, S, by = delta), y = S, lambda = lambda)
    prev[x >= (y-S) & x <= (y)] <- prev_x
  }
  prev[x > y] <- NA
  
  return(list(x = x, y = prev))
}


# that is, get the marginal prevalence by age x assuming we specify 
# lambda (single value or vector thereof)
# alpha (single value or vector thereof)
# lx (standard lifetable lx, but inside we'll get quantiles from it)
# the use the above prev_lambda_alpha_y() to get the prev function for each lifespan

prev_x_from_lambda_alpha <- function( lives_lx, 
                                      x = seq(0, omega, by = delta), 
                                      lambdas = .5, 
                                      alphas = .5, 
                                      Ss = 20, 
                                      delta = .05, 
                                      omega = 110){
  
  # vectorize to set flexible shapes by age
  # lambdas = rep(lambda, length.out = length(lives_lx))
  # alphas = rep(alpha, length.out = length(lives_lx))
  # Ss = rep(S, length.out = length(lives_lx))
  # 
  # get prevalence trajectory for each lifespan
  Prev_lifespan = matrix(0, ncol=length(x), nrow=length(lives_lx))
  for (i in 1:length(lives_lx)){
    # i = 100
    Prev_lifespan[i,] = prev_lambda_alpha_y(
                                    x = x, 
                                    y = lives_lx[i], 
                                    lambda = lambdas[i], 
                                    alpha = alphas[i], 
                                    S = Ss[i],
                                    delta = delta,
                                    omega = omega)$y
  }
  Prev_x <- colMeans(Prev_lifespan, na.rm = TRUE)
  # get mean
  #  lives = sapply(x, function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
  #  Prev_x = colSums(Prev_lifespan)/lives
  return(list(Prev_x = Prev_x, Prev_lifespan = Prev_lifespan, x = x))
}


# a function to bin prevalence into single ages
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
      prev_xout[ls, age+1] = sum(prev_x[ls, xout_index==age], na.rm = T) * delta_implied / min(1,age_d-age, na.rm = T)
    }
  }
  prev_xout
}

# funct to sum cross ages
prev_cross_x = function(Prev_lifespan_bin, lives_lx, xout){
  alives = sapply(xout, function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
  Prev_x = colSums(Prev_lifespan_bin, na.rm = T)/alives
  Prev_x
}

# flies
lambda_fun = function(x=0:100, a=.5, b=0, c=0, fly = 'poly'){
  lambdas = ifelse(fly=='exp', a * exp(b*x), 
                   x^2 * c + x * b + a)
  lambdas
}
alpha_fun = function(x=0:100, a=2, b=0, c=0, fly = 'poly'){
  alphas = ifelse(fly=='exp', a * exp(b*x), x^2 * c + x * b + a)
  alphas
}





