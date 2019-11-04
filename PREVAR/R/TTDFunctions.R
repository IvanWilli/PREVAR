# Prev function
# x decimal evaluated age
# y decimal death age
# lambda is compression parameter | higher = more compressed
#prev_sv = function(x, y, lambda=.3){ # x vector of ages at death
#  p = ifelse(is.infinite(x*exp(x*lambda)), 0, x*exp(x*lambda))/(y*exp(y*lambda))
#  p = ifelse(p>1,0,p)
#  p
#}

# more robust to higher lambda:
prev_sv = function(x, y, lambda=.3){ # x vector of ages at death
  p        <- (x * exp(lambda * (x - y)))/y
  p[p > 1 | is.nan(p)] <- 0
  p
}

# analytic integral: might speed up optimizations:
prev_sv_int <- function(y,from=0,to=y,lambda=.3){
  
  # need logical to account for lambda == 0
  if (lambda == 0){
    return(y/2 - (from^2)/(2*y) - (y-to)*(to/y) - (y-to)*(1-to/y)/2)
  }
  
  (exp(to * lambda) * 
      (-1 + to  * lambda) + 
      exp(from * lambda) * 
      (1 - from * lambda)
  ) / (exp(y * lambda) * y * lambda^2)
}

# get the ages to which specified quantiles refer.
life_bins <- function(lx, x, probs =  seq(0,1,by=.05)){
  lx <- lx / lx[1]
  splinefun(x~lx, method = "monoH.FC")(probs)
}


# a utility to help optimize for t_star, which scales the prevalence area by 
# shifting the prevalence curve (implied by alpha and lambda) left or right out of the window S
# IW: tried with analytic solution but no exist:
# https://www.wolframalpha.com/input/?i=%5Cint_0%5Et%7Bx%2Ft*%5Cexp((x-t)*%5Clambda)%7D
# Area_hat = (lambda * t_star + exp(lambda * -t_star) - 1) / (lambda^2) / t_star

get_shift = function(t_star, alpha, lambda, S){   
  # the formula depends on whether we shift left or right
  if (sign(t_star) == 1){
    Area_lambda_S <- prev_sv_int(from = t_star, to = S, y = S, lambda = lambda) + t_star
    #Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S, lambda = lambda)$value + t_star
    
  } else {
    if (sign(t_star) == -1){
      Area_lambda_S <-  prev_sv_int(from = 0, to = -t_star, y = S, lambda = lambda)
      #Area_lambda_S = integrate(f = prev_sv, lower =  (t_star), upper = S, y = S + t_star, lambda = lambda)$value
    } 
  }
  
  ((alpha * (S/2)) - (Area_lambda_S)) ^ 2
}



prev_lambda_alpha_y <- function(
  x = seq(0,omega,by=delta), 
  y = max(x), 
  lambda = 1, 
  alpha = 1, 
  S = 20, 
  omega = 110,
  delta = .05){
  
  xint <- round(x * (1/delta))
  # start with....
  # 0      y-S      y    t*
  # |______|________|______|
  #            S      t*-S
  t_star <- optimize(get_shift, interval = c(-S, S), alpha = alpha, lambda = lambda, S = S, tol = 1e-12)$minimum
  t_star <- round(t_star * 1/delta) * delta
  
  prev   <- x * 0
  #names(prev) <- x
  # TR: there is index sloppiness in here still.
  # TR: try indexing differently. Use integers?
  # yint =1200; lbint= 400; ubint= 528; t_star= 33.55; S= 40; lambda= 0.1
  yint <- round(y * (1 / delta))
  
  if (sign(t_star) == 1){
    x.i       <- seq(t_star, S, by = delta)
    x.int     <- round(rev(y - x.i) * (1 / delta))
    lbint     <- min(x.int)
    ubint     <- max(x.int)
    
    prev_x    <- prev_sv(x = x.i, y = S, lambda = lambda) 
    #	  if (sum(xint >= lbint & xint <= ubint) != sum(x.int >= 0)){
    #		  stop(paste("yint=",yint,";lbint=",lbint,";ubint=",ubint,";y=",y,";t_star=",t_star,";S=",S,";lambda=",lambda))
    #	  }
    if (t_star < y){
      prev[xint >= lbint & xint <= ubint]  <- prev_x[x.int >= 0]
    }
    prev[xint >= ubint & xint <= yint]        <- 1  
    
  }
  
  #  prev_x     <- prev_sv(x = seq(0, S + t_star, by = delta), y = S, lambda = lambda)
  #  prev[x >= (y-S-t_star) & x <= (y)]     <- prev_x
  if (sign(t_star) == -1){
    x.i        <- seq(0, -t_star, by = delta)
    x.int      <- round(rev(y - x.i) * (1 / delta))
    lbint      <- min(x.int)
    ubint      <- max(x.int)
    
    prev_x     <- prev_sv(x = x.i, y = S, lambda = lambda)
    #	  if (sum(xint >= lbint & xint <= ubint) != sum(x.int >= 0)){
    #		  stop(paste("yint=",yint,";lbint=",lbint,";ubint=",ubint,";y=",y,";t_star=",t_star,";S=",S,";lambda=",lambda))
    #	  }
    prev[xint >= lbint & xint <= ubint ]     <- prev_x[x.int  >= 0]
  }
  
  
  if (sign(t_star) == 0){
    lbint      <- round((y - S) * (1 / delta))
    prev_x     <- prev_sv(x = seq(0, S, by = delta), y = S, lambda = lambda)
    prev[xint >= lbint & xint <= yint]            <- prev_x
  }
  prev[xint > yint] <- NA
  
  return(list(x = x, y = prev))
}

#plot(prev_lambda_alpha_y(seq(0,50,by=delta),lambda=0,y=50,alpha=.1),type='l',ylim=c(0,1))
# TR: this was mostly problem withh logical matching of digits that were
# obtained via different operations. Straight 'R inferno'. 
# IW -> to get the warnings: print combinations with error, but still no finding the pattern
# for (lambda in seq(0,2,.1)){
#   for (alpha in seq(0,2,.1)){
#     for (S in seq(5, 20, 5)){
#       # lambda=1; alpha=1.5; S=10
#       t_star <- optimize(get_shift, interval = c(-S, S), alpha = alpha, lambda = lambda, S = S, tol = 1e-12)$minimum
#       t_star <- round(t_star * 1/delta) * delta
#       prev   <- x * 0
#       if (sign(t_star) == 1){
#         prev_x    <- prev_sv(x = seq(t_star, S, by = delta), y = S, lambda = lambda) 
#         x_implied <- seq((y - S), (y - t_star), by = delta)
#         tryCatch(prev[x >= (y - S) & x <= (y - t_star)] <- prev_x[x_implied >= 0], 
#                  warning=function(w) print(c(alpha,S,lambda)))
#       }}}}


# that is, get the marginal prevalence by age x assuming we specify 
# lives_lx the ages to which lifetable quantiles refer
# x in delta intervals
# lambdas (single value or vector thereof)- controls compression
# alphas (single value or vector thereof) - area relative to S (1 = S/2)
# Ss reference window width
# the use the above prev_lambda_alpha_y() to get the prev function for each lifespan

prev_x_from_lambda_alpha <- function( lives_lx, 
                                      x = seq(0, omega, by = delta), 
                                      lambdas = .5, 
                                      alphas = .5, 
                                      Ss = 20, 
                                      delta = .05, 
                                      omega = 110){
  n <- length(lives_lx)
  if (length(lambdas) == 1){
    lambdas <- rep(lambdas, n)
  }
  if (length(alphas) == 1){
    alphas <- rep(alphas, n)
  }
  if (length(Ss) == 1){
    Ss <- rep(Ss, n)
  }
  
  
  # get prevalence trajectory for each lifespan
  Prev_lifespan <- matrix(0, ncol = length(x), nrow = length(lives_lx))
  for (i in 1:length(lives_lx)){
    # i = 100
    Prev_lifespan[i,] <- prev_lambda_alpha_y(
      
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

# TR: OK if it works: this is way more complicated than I was imagining
# what about:

bin_prev_x_tr <- function(prev_x,xin){
  xout <- floor(xin)
  tapply(prev_x, xout, mean, na.rm = TRUE)
}
# IW: seems needs some fix this simplest way. Also I dont know if takes account age deaths not integers
# tried both:
# prev_x = prev_x_from_lambda_alpha(lives_lx, x = seq(0, omega, by = delta), lambdas = .5, alphas = .5, Ss = 20, delta = .05, omega = 110)$Prev_lifespan
# bin_prev_x_tr(pp, xin = x)
# bin_prev_x(pp)


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


# TR: there's a function colMeans() that would make this easier. No need for denom, right?

# IW: I think doesn't take care of (lives_lx-x) time-persons. Maybe I'm wrong

# funct to sum cross ages
prev_cross_x = function(Prev_lifespan_bin, lives_lx, xout){
  alives = sapply(xout, function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
  Prev_x = colSums(Prev_lifespan_bin, na.rm = T)/alives
  Prev_x
}

# function that do all the process to get a prev by age (useful to optimize)
prev_x_final = function(lives_lx, x, lambdas, alphas, Ss, delta, omega, xout){
  
  # get prev matrix
  Prev_lifespan = prev_x_from_lambda_alpha(lives_lx, x, lambdas, alphas, Ss, delta, omega)$Prev_lifespan
  
  # group by integer ages
  Prev_lifespan_bin = bin_prev_x(Prev_lifespan, xin = x, xout = xout)
  
  # person-year survivors in each integer age (decimal y gives the portion also)
  Prev_x = prev_cross_x(Prev_lifespan_bin, lives_lx, xout)
  
  Prev_x
}


# families
lambda_fun = function(a=.5, b=0, family = 'lineal', x=0:100){
  if(family=='exp') {
    lambas = a * exp(b*x)
  }
  if(family=='lineal'){
    lambas = x * b + a
  } 
  lambas
}

alpha_fun = function(a=2, b=0, family = 'lineal', x=0:100){
  if(family=='exp') {
    alphas = pmin(a * exp(b*x), 2)
  }
  if(family=='lineal'){
    alphas = pmin(x * b + a, 2)
  } 
  alphas
}

#plot(alpha_fun(.5,.01, family="lineal"))
##  1024 64  0 4 64 144 64

expit <- function(x){
  exp(x) / (1 + exp(x))
}


# started experimenting with polynomials
my_poly <- function(pars,x){
  n        <- length(pars)
  orders   <- 1:n - 1
  rowSums(outer(x,orders,"^") %*% diag(pars) )
}

# alpha ranges from 0-2, so 2*expit should guarantee that.
alpha_poly <- function(x,...){
  pars <- unlist(list(...))
  out  <- 2 * expit(my_poly(pars, x))
  out[is.nan(out)] <- 2
  out[is.na(out)] <- 0
  out
}

lambda_poly <- function(x,...){
  pars <- unlist(list(...))
  # strictly positive
  out <- exp(my_poly(pars, x))
  out[is.nan(out)] <- 0
  out[is.na(out)] <- 0
  out
}

# x  <- seq(.5,1,length=100)
# p1 <- seq(0,10,length=10)
# p2 <- seq(0,10,length=10)
# 
# plot(x,dbeta(x,12,5))



beta_min <- function(pars = c(scale = 1, shape1 = 10, shape2 = 5), 
                     lx, x, piObs){
  
  xbeta   <- seq(.5,1,length = length(x)) # ^ pars["z"]
  
  lambdas <- pars["scale"] * dbeta(xbeta, shape1 =  pars["shape1"], shape2 =  pars["shape2"])
  
  PiEst   <- getPiEst(lx,x,lambdas)
  
  sum((PiEst - piObs)^2 )
}




