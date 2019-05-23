# Survival from HMD data
library(HMDHFDplus)
us = 'ivanwilliams1985@gmail.com'; pw = 'volveroman'
mlt <- readHMDweb("USA","mltper_1x1", us, pw)
lx <- subset(mlt,Year == 2010, "lx")$lx

# get ages when survival is 1, .99, .98, ..., .01, 0
q <- seq(1,0,by=-.01)
lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
plot(lives_lx, q[-1] , type = 's')

# show some different parameter combinations
plot(NULL, type = "n", xlim = c(0,50), ylim = c(0,1))
lines(prev_lambda_alpha_y(y=50,alpha=1,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=.5,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=1.5,lambda=0,S=20))
lines(prev_lambda_alpha_y(y=50,alpha=.5,lambda=.5,S=20))

# get matrix LSxAges of prevalences
Prev_lifespan <- prev_x_from_lambda_alpha(lx)$Prev_lifespan

# For fixed rewards the variance is *zero* unless alpha, lambda, S vary over age.
rowSums(Prev_lifespan,na.rm=TRUE) * delta
# for binomal variance we'd need to invert the prevalence curve 
# in each lifespan. Make sense? Or at least do a discrete inversion,
# in case a prevalence doesn't make it to 1 by death.
# that is a to-do.

# -------------------------------------------------
# Now make a demonstrative plot for IW:
omega <- 110
q        <-  seq(1, 0, by = -.05)
lives_lx  <- round(life_bins(lx, 0:110, probs = q)[-1])
lives_lx[lives_lx>omega] <- omega
xx        <- seq(0, 110, by = .1)
S         <- lives_lx * .1 # 10% of life assume
alphas    <- seq(.1, .5, length = length(lives_lx))
lambdas   <- seq(5, .5, length = length(lives_lx))


PrevMat <- matrix(ncol = 1101,nrow=20)

plot(NULL, type = 'n', xlab = "Age", ylab = "Proportion", axes = FALSE, xlim = c(0,110), ylim=c(0,1))
rect(0,q[-length(q)],lives_lx,q[-1], border = gray(.8))
lines(c(0,rep(lives_lx,each=2),0), rep(q,each=2))

for (i in 1:length(lives_lx)){
  pxi <- prev_lambda_alpha_y(
    x = xx, 
    y = lives_lx[i], 
    lambda = lambdas[i], 
    alpha = alphas[i], 
    S = S[i],
    delta = .1,
    omega = omega)
  xx <- pxi$x
  px <- pxi$y
  # scale down
  px <- px * .05
  PrevMat[i, ] <- px
  rem <- is.na(px) | px == 0
  ends <- range(xx[!rem])
  polygon(
    x=c(ends[1],xx[!rem],ends[2]),
    y=c(q[i+1], px[!rem] + q[i+1], q[i+1]) ,
    col = "black"
  )
}

prev_L <- rowSums(PrevMat, na.rm=TRUE) 
var_L = sum((prev_L - mean(prev_L)) ^ 2) / length(prev_L)
var_L^.5/mean(prev_L)
# ------------------------

Prev_lifespan_bin = bin_prev_x(prev_x_from_lambda_alpha(lx)$Prev_lifespan, 
                               xin = seq(0,110,by=.05), 
                               xout = 0:100)
