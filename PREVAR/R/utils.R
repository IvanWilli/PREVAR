
# Author: tim
###############################################################################
qx2lx <- function(qx,radix=1e5){
	n <- length(qx)
	px <- 1-qx
	radix * c(1,cumprod(px))[1:n]
}

lx2qx <- function(lx){
	lx <- c(lx, 0)
	lx <- lx / lx[1]
	n  <- length(lx)
	px <- lx[-1] / lx[-n]
	1-px
}

lx2dx <- function(lx,radix=1e5){
	lx <- lx / lx[1]
	-diff(radix * c(lx,0))
}

granularize <- function(y_single, 
		int_x = 1:length(y_single) - 1, 
		interval = .5, 
		method = c("monoH.FC","linear")[1],
		xsqrt = TRUE){
	
	xnew <- seq(min(int_x),max(int_x),by=interval)
	if (xsqrt){
		xnew  <- sqrt(xnew)
		int_x <- sqrt(int_x)
	}
	if (method == "linear"){
		out <- approx(x=int_x, y=y_single, xout = xnew)$y
	}
	if (method == "monoH.FC"){
		out <- splinefun(x = int_x, y = y_single, method = "monoH.FC")(xnew)
	}
	out
}

prev_line <- function(vec, from = 0, to = .5){
	seq(from,to,length=length(vec))
}

'%==%' <- function(x,y){
	if (is.null(x)){
		return(rep(FALSE, length(y)))
	}
	if (is.null(y)){
		return(rep(FALSE,length(x)))
	}
	x == y & !is.na(x) & !is.na(y)
}


'%!=%' <- function(x,y){
	if (is.null(x)){
		return(rep(TRUE, length(y)))
	}
	if (is.null(y)){
		return(rep(TRUE,length(x)))
	}
	x != y & !is.na(x) & !is.na(y) 
}