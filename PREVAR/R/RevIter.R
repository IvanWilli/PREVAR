lives_lx = c(0.8, 1.1, 2.8, 3.9, 4.7, 5.1, 6.05, 7.6, 8.4, 10)
x = 0:10
S = pmin(3,lives_lx)
Pi_Obs = round(.8/(1+exp(-.1*(x[-length(x)]-7))),2)
barplot(rev(lives_lx), horiz = T, col = 0, space = 0, xlim = c(0,max(x)), names.arg = rev(lives_lx))
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
    # i = 4
    lives_i = groups[i,][!is.na(groups[i,])]
    ages = rounds[i,][!is.na(rounds[i,])]
    
    lambda[lives_i] = optimize(get_lambda, lower = -50, upper = 50, 
                               x_ini = min(ages), x_fin = lives_lx[lives_i],
                               PrevObs_x = sum(Pi_Obs[ages+1] * alives[ages+1]) - 
                                           ifelse(i==nrow(groups), 0, sum(prevs[max(lives_i):nrow(prevs), ages+1])), 
                               S = S[lives_i], 
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
Pi_Obs * alives; colSums(prevs)
round(Pi_Obs * alives[-6] - colSums(prevs),5)
sum((Pi_Obs * alives)[3:5])
sum(colSums(prevs)[3:5])

# graph
for(s in 1:length(lives_lx)){
  lines(seq(0,lives_lx[s],.05), prev_sv_simetric(seq(0,lives_lx[s],.05), 
                                                 y = lives_lx[s], S = S[s], 
                                                 lambda = lambda[s])+(length(lives_lx)-s))
}
plot(lambda)

# get lambda function
get_lambda = function(lambda, S, x_ini, x_fin, PrevObs_x, lives){
  Prev_x = c()
  for(j in 1:lives){
    Prev_x[j] = integrate(f = prev_sv_simetric, lower =  x_ini, upper = x_fin[j], subdivisions = 20000, 
                       y = x_fin[j], S = S[j], lambda = lambda)$value
  }
  (PrevObs_x - sum(Prev_x))^2
}

















