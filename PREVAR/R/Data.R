
# Author: tim
###############################################################################

install.packages("HMDHFDplus")
library(HMDHFDplus)
LT        <- readHMDweb("USA","mltper_1x1",username=us,password=pw)
#lx   <- LT$lx[LT$Year==2000]
qx        <- LT$qx[LT$Year==2000]
