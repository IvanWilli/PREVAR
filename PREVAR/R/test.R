############################ OPtim with constraints
# Ages >= lambdas: mas ecuaciones q incognitas

# Prev function
# x decimal evaluated age
# y decimal death age
# lambda is compression parameter | higher = more compressed
prev_sv = function(x, y, lambda=.3){ # x vector of ages at death
	p = x*exp(x*lambda)/(y*exp(y*lambda))
	p = ifelse(p>1,0,p)
	p
}

# invent data
y = c(1, 5, 8)
n = m = length(y) # same ages at death than individuals
Ages = 0:max(y-1)
Pi_obs = data.frame(x = Ages, Pi = .0005*exp(Ages*.8)) # at begining of ages
lives = sapply(Ages, function(x) sum(y>x))
barplot(Pi_obs$Pi, names.arg = Ages, space = 0, xlab = 'Age', ylab = 'Prev')
text(rep(.1,length(Ages)),labels=as.character(lives), adj = c(6,1), col=2)

# obj funct
fun_opt2 = function(lambda){
	Pi_hat_Age = c()
	for (age in Ages){
		# age=3
		Pi_hat_Age_i = sum(sapply(1:n, function(i) integrate(f = prev_sv, lower =  age, upper = age+1,
									y = y[i], lambda = lambda[i])$value))
		Pi_hat_Age = c(Pi_hat_Age, Pi_hat_Age_i)
	}
	sum((Pi_obs$Pi - Pi_hat_Age/lives)^2)
}

# ineq
fun_ineq = function(lambda){
	z1 = lambda[1]-lambda[2]
	z2 = lambda[2]-lambda[3]
	return(c(z1,z2))
}

# iter optim
#install.packages("Rsolnp")
library(Rsolnp)
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,max(y)), xlab='Age', ylab='prev')
points(Pi_obs$x+.5, Pi_obs$Pi, col=4, lty=2, lwd=2, t='o')
lambda_hats = Prev_hat = matrix(0,nrow=1,ncol=length(Ages))
for (s in 1:20){
	lambda0 = sample(seq(1,8,.5),n)
	# optim
	optim = solnp(lambda0, fun = fun_opt2,ineqfun = fun_ineq,
			ineqLB = c(0,0), ineqUB = c(10,10), LB = c(0,0,0), UB = c(50,50,50),
			control = list(delta = .1))
	lambda_hat = optim$pars
	lambda_hats = rbind(lambda_hats, lambda_hat)
	x = seq(0, max(y), .05)
	for (i in 1:length(lambda_hat)){
		lines(x[x<=y[i]], prev_sv(x[x<=y[i]], y[i], lambda_hat[i]),col=i)
	}
	Pi_hat_Age = c()
	for (age in Ages){
		Pi_hat_Age_i = sum(sapply(1:n, function(i) integrate(f = prev_sv, lower =  age, upper = age+1,
									y = y[i], lambda = lambda_hat[i])$value))
		Pi_hat_Age = c(Pi_hat_Age, Pi_hat_Age_i)
	}
	Prev_hat =  rbind(Prev_hat, Pi_hat_Age/lives)

}

Prev_hat
lambda_hats = lambda_hats[-1,]
# error
Pi_obs_matrix = matrix(rep(Pi_obs$Pi,s), nrow=s, ncol=length(Ages), byrow=T)
error.Abs_iter = Prev_hat[-1,] - Pi_obs_matrix
error.Relat_iter = round(error.Abs_iter/Pi_obs_matrix*100,2)
round(Pi_obs$Pi,4)
round(colMeans(Prev_hat[-1,]),4)
colMeans(abs(error.Abs_iter))
colMeans(abs(error.Relat_iter))

Prev_hat

life_bins <- function(lx, x, probs =  seq(0,1,by=.05)){
	lx <- lx / lx[1]
	splinefun(x~lx, method = "monoH.FC")(probs)
}

library(HMDHFDplus)
mlt <- readHMDweb("USA","mltper_1x1",us,pw)
lx <- subset(mlt,Year == 2010, "lx")$lx


lives_lx <- life_bins(lx, 0:110, probs = seq(1,0,by=-.01))

plot(lives, seq(1,0,by=-.1), type = 's' )

delta <- .05
round(lives_lx / delta) * delta

l2 <- 17.29474
x <- seq(0,2,by=.05)
plot(x,prev_sv(x, 2, l2), col = "red", type = 'l')
lines(splinefun(x~prev_sv(x, 2, l2), method = "monoH.FC")(seq(1,0,by=-.01)),seq(1,0,by=-.01), col = "blue", lty = 2)


dis_x <- 2-splinefun(x~prev_sv(x, 2, l2), method = "monoH.FC")(seq(1,.01,by=-.01))
sum((dis_x - mean(dis_x)) ^ 2)

