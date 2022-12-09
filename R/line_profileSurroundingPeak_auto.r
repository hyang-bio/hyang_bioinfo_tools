library(RColorBrewer)

args <- commandArgs(T)
inF <- args[1]
outFL <- args[2]
widthI <- as.double(args[3])
heightI <- as.double(args[4])
xL <- strsplit(args[5], split = ';')[[1]] # -30kb;HMR_S;HMR_T;30kb
xLoc <- as.numeric(strsplit(args[6], split = ';')[[1]]) # 1;30;45;75
yL <- args[7]
mainL <- args[8]
colN <- args[9]


data <- read.table(inF, sep = '\t', header = F)
if(grepl(';', colN)){
	colorRamp <- strsplit(colN, split = ';')[[1]]
	}else{
		colorRamp <- colorRampPalette(brewer.pal(n = 7, name = colN))(nrow(data))
	}

data[, 2:ncol(data)] <- abs(data[, 2:ncol(data)])

ylim_l <- min(data[, 2:ncol(data)])-0.15*abs(min(data[, 2:ncol(data)]));print(ylim_l)
ylim_u <- quantile(unlist(data[, 2:ncol(data)]), 0.9995)+0.15*abs(quantile(unlist(data[, 2:ncol(data)]), 0.9995));print(ylim_u)

pdf(paste0(outFL, '.pdf'), width = widthI, height = heightI)
plot(seq(1, ncol(data)-1, 1), data[1, 2:ncol(data)], xlim = c(0, ncol(data)-1), ylim = c(ylim_l, ylim_u), xlab = '', ylab = yL, xaxt = 'n', main = list(mainL, cex = 0.9), type = 'l', lwd = 1, col = colorRamp[1])
if(nrow(data) >= 2){
	for(i in 2:nrow(data)){
		lines(seq(1, ncol(data)-1, 1), data[i, 2:ncol(data)], lwd = 1, lty = 1, col = colorRamp[i])
	}
}
axis(side = 1, xLoc, labels = xL)
legend('topright', legend = data[, 1], col = colorRamp, lty = 1, lwd = 2, bty = 'n', cex = 0.5)
dev.off()