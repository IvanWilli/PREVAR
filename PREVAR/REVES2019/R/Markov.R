
# Author: tim
###############################################################################

setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
source("R/BlockFunctions.R")
#----------------------

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx), interval = .05, method = "monoH.FC"))  
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05, 0, .5)

var_fixed <- VarFixed(qx_05,prev_05,interval=.05,calcs = "faster",closeout=TRUE)
var_bern  <- VarBernoulli(qx_05,prev_05,interval=.05,calcs = "faster",closeout=TRUE)

saveRDS(c(bern=var_bern,fixed=var_fixed),"Data/varbernfixed.rds")

###########################################################
# need code to draw fixed rewards prevalence, same 5-year 100 radix setup.