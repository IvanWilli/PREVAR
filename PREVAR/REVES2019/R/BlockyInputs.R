
# Author: tim
###############################################################################
setwd("/home/tim/git/PREVAR/PREVAR")
source("R/Data.R")
source("R/utils.R")
source("R/PrevRewards.R")
source("R/BlockFunctions.R")
source("R/BlockFunctions.R")

lx         <- qx2lx(qx)
lx_05      <- exp(granularize(log(lx),interval=5,method = "monoH.FC"))
qx_05      <- lx2qx(lx_05)
prev_05    <- prev_line(lx_05,0,.5)
prev       <- prev_line(lx,0,.5)
a5         <- seq(0,110,by=5)
a          <- 0:110


# original lifetable and prevalence:

pdf("REVES2019/Figs/Inputs.pdf")
par(xaxs = "i", yaxs = "i")
plot(a, lx / lx[1], type = 'l', las = 1,
		xlab = "Age",
		ylab = "Probability",
		cex.lab = 1.5,
		cex.axis = 1.5,
		lwd = 3,
		ylim = c(0, 1.01),
		panel.first=list(
				polygon(x = a, y = prev * lx / lx[1], col = "#FF000080", border = NA)))
lines(a, prev, col = "red", lwd = 3)
dev.off()

pdf("REVES2019/Figs/InputsBlocky.pdf")
par(xaxs = "i", yaxs = "i")
plot(a5, lx_05 / lx_05[1] * 100, type = 'n', las = 1,
		xlab = "Age",
		ylab = "100 lives",
		cex.lab = 1.5,
		cex.axis = 1.5,
		lwd = 3,
		ylim = c(0, 101),
		panel.first = list(
				polygon(c(0,rep(a5[-1],each=2),110),
						round(rep(prev_05 * 100 * lx_05 / lx_05[1],each=2)),
						border = NA,
						col = "#FF000080")))
lx_block_line(lx_05, a5, radix = 100, lwd = 3, round = TRUE)
lines(c(0,rep(a5[-1],each=2),110), rep(prev_05*100,each=2), col = "red", lwd = 3)
dev.off()

