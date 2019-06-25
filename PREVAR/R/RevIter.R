get_lambda_sim = function(lambda, S, x_ini, x_fin, PrevObs_x, lives){
	# x_ini=8; x_fin=c(8.7); S=5; lambda=.5; lives=1
	Prev_x = c()
	for(j in 1:lives){
		Prev_x[j] = integrate(f = prev_sv_simetric, lower =  x_ini, upper = x_fin[j], subdivisions = 20000, 
				y = x_fin[j], S = S[j], lambda = lambda)$value
	}
	(PrevObs_x - sum(Prev_x))^2
}

# prev function
prev_sv_simetric = function(x, y, S, lambda=3){ # x vector of ages at death
	# x = c(1.5, 2, 2.5, 3, 3.5, 3.7, 4, 4.5); y = 4; S =2; lambda = -1
	xS = pmax(x-(y-S), 0)
	p <- (xS * exp(abs(lambda) * (xS - S)))/S
	p[p > 1 | is.nan(p)] <- 0
	if(lambda<0){
		xS = pmax(S-(x-(y-S)), 0)
		p <- (xS * exp(abs(lambda) * (xS - S)))/S
		p[x>=(y-S) & x<=y] = 1-p[x>=(y-S) & x<=y]
		p[p > 1 | is.nan(p)] <- 0
	}
	p
}
plot(0,0, col=0, ylim = c(0,1), xlim=c(0,10), xlab='x', ylab='prev', asp=T)
for (lambdas in seq(-10, 10, .5)){
	lines(seq(0, 10, .05), prev_sv_simetric(seq(0, 10, .05), y = 10, S = 5, lambdas), 
			col=ifelse(lambdas>0, 1, 2))
}

# ok
#lives_lx = c(0.8, 1.6, 2.9, 3.8, 4.7, 6.3, 6.56, 7.6, 8.7, 10)
# no ok
 lives_lx = c(0.8, 1.6, 2.1, 3.8, 4.7, 6.3, 6.56, 7.6, 8.7, 10)
x = 0:10
# q = seq(1, 0, -.1)
# lives_lx <- life_bins(lx, 0:110, probs = q )[-1]
# lives_lx[lives_lx>=100] = 100
# lives_lx = lives_lx[lives_lx>50]
# x = 50:100
S = pmin(5,lives_lx-min(x))

Pi_Obs = .03*exp((x[-length(x)]+.5)*.3)
# Pi_Obs = .03*exp((x[-length(x)]+.05)*.03)
plot(x[-length(x)], Pi_Obs, t='o')
barplot(rev(lives_lx)-min(x), horiz = T, col = 0, space = 0, xlim = range(x)-min(x), names.arg = rev(lives_lx))
abline(v=1:max(x), lty=2)
mtext(text = round(Pi_Obs,2), side = 1, line = 2, at = seq(.6, max(x)-.4, 1), col = 2)
mtext('Prev', side = 1, line = 2, at = -.5, col = 2, cex=.8)

# how many PERSON-YEARS at each x
alives = sapply(x[-length(x)], function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
mtext(text = alives, side = 1, line = 3, at = seq(.6, max(x)-.4, 1), col = 4)
mtext('Pers-years', side = 1, line = 3, at = -.5, col = 4, cex=.8)

# groups and rounds

groups0 = unique(sapply(x[-length(x)], FUN = function(s) which(lives_lx - s>0)))
groups = list()
ncol_groups = 1
for(s in length(groups0):1){
  groups[[s]] = groups0[[s]]
  if(s < length(groups0)){
    groups[[s]] = setdiff(groups0[[s]],groups0[[s+1]])
  }
  if(length(groups[[s]])==0){
    groups[[s]] = groups[[s+1]]
  }
  ncol_groups = max(ncol_groups, length(groups[[s]]))
}
groups_mat = matrix(NA,length(groups), ncol_groups)
for(s in 1:length(groups)){
  groups_mat[s,1:length(groups[[s]])] = groups[[s]]
}
groups = groups_mat

rounds0 = sapply(x[-length(x)], FUN = function(s) lives_lx - s>0)
rounds0 = unique(sapply(1:nrow(rounds0), function(x) which(rounds0[x,]==T)))
rounds = list()
ncol_rounds = 1
for(s in 1:length(rounds0)){
  rounds[[s]] = rounds0[[s]]
  if(s > 1){
    rounds[[s]] = setdiff(rounds0[[s]], rounds0[[s-1]])
  }
  if(length(rounds[[s]])==0){
    rounds[[s]] = rounds[[s-1]]
  }
  ncol_rounds = max(ncol_rounds, length(rounds[[s]]))
}
rounds_mat = matrix(NA,length(rounds), ncol_rounds)
for(s in 1:length(rounds)){
  rounds_mat[s,1:length(rounds[[s]])] = rounds[[s]]-1
}
rounds = rounds_mat

# lambdas y matriz de prevalencias
prevs = matrix(0, length(lives_lx), length(x)-1)
colnames(prevs) = x[-length(x)]
rownames(prevs) = round(lives_lx,2)
lambda = c()

# which lambda replicates prev in tha
for(i in nrow(groups):1){
    # i = 8
    lives_i = groups[i,][!is.na(groups[i,])]
    ages = rounds[i,][!is.na(rounds[i,])]
    
    
    
    
    lambda[lives_i] = optimize(get_lambda_sim, lower = -10, upper = 10, 
                               x_ini = x[min(ages)+1], x_fin = lives_lx[lives_i],
                               PrevObs_x = sum(Pi_Obs[ages+1] * alives[ages+1]) - 
                                           ifelse(i==nrow(groups), 0, 
                                                  sum(prevs[max(lives_i):nrow(prevs), ages+1])), 
                               S = S[lives_i], 
                               lives = length(lives_i), tol = 1e-10)$minimum
    
    for(live in lives_i){
      prevs[live,] = sapply(x[-length(x)], 
                            FUN = function(s) integrate(f = prev_sv_simetric, lower =  s, upper = s+1,  
                                                                       y = lives_lx[live], S = S[live], 
                                                                       lambda = lambda[live],
                                                                       subdivisions = 20000)$value)
      
    }
}

# graph
for(s in 1:length(lives_lx)){
  lines(seq(0,lives_lx[s],.05), prev_sv_simetric(seq(0,lives_lx[s],.05), 
                                                 y = lives_lx[s], S = S[s], 
                                                 lambda = lambda[s])+(length(lives_lx)-s))
}

#chek!
fit = round(matrix(c(Pi_Obs * alives, colSums(prevs)), nrow = 2, ncol = length(alives), byrow = T), 2)
barplot(fit, beside = T, names.arg = x[-length(x)], legend.text = c('Obs', 'Fit'), ylab = 'PY prevalence')
plot(lambda)

# get lambda sim function

