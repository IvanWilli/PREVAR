
# Author: tim
###############################################################################


makeLTblock <- function(lx, prev, radix, interval){
	n       <- length(lx)
	stopifnot(length(prev) == n)
	lx      <- lx / lx[1] * radix
	Nsick   <- round(lx * prev)
	LTBlock <- matrix(0, ncol = length(lx), nrow = radix)
	for (i in 1:n){
		if (Nsick[i] > 0){
			ind <- sample(1:lx[i], size = Nsick[i], replace = FALSE)
			LTBlock[ind,i] <- interval
		}
	}
	LTBlock
}

# variance simulation
vari <- function(lx, prev, radix = 1e5, interval=1){
	LTBlock <- makeLTblock(
			lx = lx, 
			prev = prev, 
			radix = radix, 
			interval = interval)
	Di      <- rowSums(LTBlock)
	DLE     <- mean(Di) # same every time!
	mean((Di - DLE)^2)
}

