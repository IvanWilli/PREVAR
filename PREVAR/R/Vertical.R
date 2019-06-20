# x_len is the grid in age
# lx_len is the height of each vertical bar
# S is a decreasing function of x
# Fisrt estimate lambda for each vertical bar, and then transpose to have the final matrix

####### prev function
prev_sv_simetric = function(x, y, S, lambda=3){ # x vector of ages at death
  xS = pmax(x-(y-S), 0)
  p <- (xS * exp(abs(lambda) * (xS - S)))/S
  p[p > 1 | is.nan(p)] <- 0
  if(lambda<0){
    pi = sort(1 - p[p>0])
    p[p>0] = pi
  }
  p
}
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,100), xlab='x', ylab='prev')
for (lambdas in seq(-2, 2, .05)){
  lines(seq(0, 100, .05), prev_sv_simetric(seq(0, 100, .05), y = 100, S = 20, lambdas), 
        col=ifelse(lambdas>0, 1, 2))
}

####### get area prev, scaling by x_len
get_area = function(x_len, lx_len, S, lambda){
  integrate(f = prev_sv_simetric, lower =  0, upper = lx_len, y = lx_len, S = S, lambda = lambda, subdivisions = 5000)$value * x_len
}

####### get lambda for each vertical bar: x_len * lx_len
get_lambda = function(lambda, x_len, lx_len, PrevObs_x, S){
  Prev_hat_x = get_area(x_len, lx_len, S, lambda)
  (PrevObs_x * lx_len - Prev_hat_x)^2
}
# works?
lambda_optim = optimize(get_lambda, lower = -50, upper = 50, x_len = .05, lx_len = .8, S = .2, PrevObs_x = .002, tol = 1e-6)$minimum
get_area(x_len = .05, lx_len = .8, S = .2, lambda = lambda_optim) - .002 * .8

#######  S as a portion of lx_len, increasing by age
S_fun = function(x, ages = c(0,100), initial=.1, final=1) initial + (x-ages[1])/(ages[2]-ages[1]) * (final-initial)

####### Data
omega <- 110
delta <- .5
x = seq(0, omega, delta)
library(splines)




spline_lx <- interpSpline(0:omega, lx/lx[1])
lx_x <- predict(spline_lx, x)$y
plot(x, lx_x, t='s',ylab='lx',xlab='x')
Pi_Obs = function(x) .8/(1+exp(-.1*(x-70)))
lines(x, Pi_Obs(x), col = 2)
lines(x, S_fun(x, ages=c(0, omega)), col=3)
legend(5,.7,c('lx', 'prev', 'S (portion of lx)'), col=1:3, lty=1, bty = 'n')
lines(x,lx_x - lx_x * S_fun(x,ages=c(0, omega)))
####### get lambda for each delta x (vertical bar)
lambdas = c()
for(i in 1:length(x)){
  # i = 1
  lambdas[i] = optimize(get_lambda, lower = -500, upper = 500, x_len = delta, lx_len = lx_x[i], 
                        S = S_fun(x[i], ages=c(0, omega)) * lx_x[i],
                        PrevObs_x = Pi_Obs(x[i] + delta/2), 
                        tol = 1e-6)$minimum
}
plot(lambdas)
# see a case
i = 180
x[i]; lx_x[i]; S_fun(x[i], ages=c(0, omega)); lambdas[i]

###### get prev for each cell (slooow for delta -> 0)
Prev_mat = matrix(NA, length(x), length(lx_x))
for (j in 1:(length(lx_x))){
  # j = 1
  Prev_mat[j,i] = integrate(f = prev_sv_simetric, lower =  lx_x[j], upper = lx_x[j+1], y = lx_x[j], S = lx_x[j]*.8, 
                            lambda = lambda_i, 
                            subdivisions = 5000)$value * x[i]
}

### transpose to get straight matrix

Prev_mat = rev(t(Prev_mat))
heatmap(Prev_mat, Rowv = NA, Colv = NA)



lx_inv <- function(lx, x, probs, deltax = 0){
	lx <- lx / lx[1]
	lxq <- splinefun(x ~ lx, method = "monoH.FC")(probs)
	if (delta > 0){
		lxq <- lxq - lxq %% deltax
	}
	lxq
}

lx_inv(lx,0:110,seq(1,0,by=-.001),deltax=.05)































lambda_optim = optimize(get_lambda, lower = -50, upper = 50, x_len = .5, lx_len = .8, PrevObs_x = .02, 
                        tol = 1e-6)$minimum
prev = prev_sv_simetric(x = seq(0,lx_len,.01), y = lx_len, S = lx_len/2, lambda = lambda_optim) * x_len
plot(seq(0, lx_len, .01), prev, t='l')

y <- c(71, 34, 11, 9.6, 26, 50, 38, 38, 30, 36, 31)
n <- length(y)
x <- 1:n
s = smooth.spline(x, y, spar=0.5)
xy <- predict(s, seq(min(x), max(x), by=1)) # Some vertices on the curve
m <- length(xy$x)                         

x.poly <- c(xy$x, xy$x[m], xy$x[1])         # Adjoin two x-coordinates
y.poly <- c(xy$y, 0, 0)                     # .. and the corresponding y-coordinates
plot(range(x), c(0, max(y)), type='n', xlab="X", ylab="Y")
polygon(x.poly, y.poly, col=gray(0.95), border=NA)          # Show the polygon fill only
lines(s)
points(x.poly, y.poly, pch=16, col=rainbow(length(x.poly))) # (Optional)