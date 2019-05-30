
# Author: tim
###############################################################################
var_rand <- readRDS("Data/vari1000_1e5_05.rds")
limits   <- readRDS("Data/blockminmax.rds")
markov   <- readRDS("Data/varbernfixed.rds")

var_comp <- sort(c(random = mean(var_rand),
		min = limits["min"], 
		max = limits["max"],
		top = limits["top"],
		bern = markov["bern"],
		fixed = markov["fixed"]))

pdf("REVES2019/Figs/Comparison.pdf")
dotchart(sqrt(var_comp),
		labels = c("Minimum","Random","Fixed","Bernoulli","Top","Max (bottom)"),
		cex = 1,
		main = "", 
		xlab = "sd",
		pch = 16)
dev.off()