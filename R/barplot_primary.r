# usage: Rscript barplot.r dataFile outFile mainLabel
# Notes: input file is separated by ','


args <- commandArgs(T)
dataFile <- args[1]
mainLabel <- args[2]
outFile <- args[3]
ylim_l <- as.double(args[4])
ylim_u <- as.double(args[5])
breaks <- as.double(args[6])
pdf_width <- as.integer(args[7])
pdf_height <- as.integer(args[8])
# print(args)


cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")

if(file.info(dataFile)$size != 0){
	data <- read.table(dataFile,sep='\t',row.names=1,header = T)
	pdf(outFile, width = pdf_width, height = pdf_height)
	bp <- barplot(as.numeric(data[,1]),col=cccol[1:dim(data)[1]],ylab='Enrichment Score',ylim = c(ylim_l,ylim_u), yaxt = 'n',main=mainLabel,las=2, names = rownames(data), cex.names = 0.75)
	abline(h = -1, lty = 2, lwd = 2, col = 'gray');abline(h = 1, lty = 2, lwd = 2, col = 'gray')
	axis(side = 2, at = seq(ylim_l,ylim_u,by=breaks))
	# text(x = bp, y = as.double(data[,1]*0.9999), labels = data[,1], pos = 1)
	box(lwd = 2)
	dev.off()
}
