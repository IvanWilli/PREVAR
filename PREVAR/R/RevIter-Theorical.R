lives_lx = c(rep(1,2), 2, 4, rep(5, 2), 7, rep(8,2), 10)
x = 0:10
S = pmin(3,lives_lx)
Pi_Obs = round(.7/(1+exp(-.6*(x[-length(x)]-7))),2)
plot(Pi_Obs)
barplot(rev(lives_lx), horiz = T, col = 0, space = 0, xlim = c(0,max(x)), names.arg = rev(lives_lx), yaxt="n")
abline(v=1:max(x), lty=2)
mtext(text = Pi_Obs, side = 1, line = 2, at = seq(.6, max(x)-.4, 1), col = 2)

# how many PERSON-YEARS at each x
alives = sapply(x[-length(x)], function(x) sum(pmin(pmax(lives_lx-x, 0), 1)))
mtext(text = alives, side = 1, line = 3, at = seq(.6, max(x)-.4, 1), col = 4)

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


# rounds = matrix(c(0,NA,NA,1,2,3,4,NA,NA), 3, 3, byrow = T)
# groups = matrix(c(1, NA, 2, 3, 4, NA), 3, 2, byrow = T)

# lambdas y matriz de prevalencias
prevs = matrix(0, length(lives_lx), length(x)-1)
colnames(prevs) = x[-length(x)]
rownames(prevs) = lives_lx
lambda = c()

# which lambda replicates prev in tha
for(i in nrow(groups):1){
    # i = 1
    # que vida/s es/son
    lives_i = groups[i,][!is.na(groups[i,])]
    # de qué edades es/son responsable/s
    ages = rounds[i,][!is.na(rounds[i,])]
    
    lambda[lives_i] = optimize(get_lambda, 
                               S = S[lives_i], 
                               lower = -50, 
                               upper = 50, 
                               x_ini = min(ages), 
                               x_fin = lives_lx[lives_i],
                               PrevObs_x = sum(Pi_Obs[ages+1] * alives[ages+1]) - 
                                           ifelse(i==nrow(groups), 0, sum(prevs[max(lives_i):nrow(prevs), ages+1])), 
                               lives = length(lives_i))$minimum

    
    for(live in lives_i){
      prevs[live,] = sapply(x[-length(x)], 
                            FUN = function(s) integrate(f = prev_sv_simetric, lower =  s, upper = s+1,  
                                                                       y = lives_lx[live], S = S[live], 
                                                                       lambda = lambda[live],
                                                                       subdivisions = 20000)$value)
    }
}


#chek!
(Pi_Obs * alives - round(colSums(prevs),2))/Pi_Obs/alives

# graph
for(s in 1:length(lives_lx)){
  lines(seq(0,lives_lx[s],.05), prev_sv_simetric(seq(0,lives_lx[s],.05), 
                                                 y = lives_lx[s], S = S[s], 
                                                 lambda = lambda[s])+(length(lives_lx)-s))
}
mtext(text = round(colSums(prevs),2), side = 1, line = 4, at = seq(.6, max(x)-.4, 1), col = 3)
mtext(text = rev(round(lambda,1)), side = 2, line = .5, at = seq(.6, length(lambda)-.4, 1), 
      col = 9, cex=.6, las = .5)
mtext(text = "lambdas", side = 2, line = 1.2, at = length(lambda)/2)

######################## get lambda function
get_lambda = function(lambda, S, x_ini, x_fin, PrevObs_x, lives){
  # lambda=1.5; S=3; x_ini=8; x_fin=10; PrevObs_x=.8; lives=1
  Prev_x = c()
  for(j in 1:lives){
    # j = 1
    Prev_x[j] = integrate(f = prev_sv_simetric, lower =  x_ini, upper = x_fin[j], subdivisions = 20000, 
                          y = x_fin[j], S = S[j], lambda = lambda)$value
  }
  (PrevObs_x - sum(Prev_x))^2
}
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




















