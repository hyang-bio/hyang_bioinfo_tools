args <- commandArgs(T)
dataFN <- args[1]
mainLabel <- args[2]
xMax <- as.double(args[3])
yMax <- as.double(args[4])
zMax <- as.double(args[5])
group_order <- strsplit(args[6], split = ';')[[1]]
pdf_width <- as.double(args[7])
pdf_height <- as.double(args[8])
code_color <- strsplit(args[9], split = ';')[[1]]


library(ggplot2)
library(scatterplot3d)
# source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

data <- read.table(dataFN, sep = '\t', header = T, check.names = F, comment.char = '')
data[, 3][data[, 3] > zMax] <- zMax
data$Group <- factor(data$Group, levels = group_order)
colors <- code_color[as.numeric(data$Group)]
pdf(paste0('3D scatterplot of ', mainLabel, '.pdf'), width = pdf_width, height = pdf_height, useDingbats = F)
	s3d <- scatterplot3d(data[, 1:3], pch = 16, color = alpha(colors, 0.4), type = 'h', xlim = c(0, xMax), zlim = c(0, zMax), main = list(mainLabel, cex = 1))
	# addgrids3d(data[, 1:3], grid = c('xy', 'xz', 'yz'))
	legend(s3d$xyz.convert(max(data[, 1])*0.8, max(data[, 2])*0.8, max(data[, 3])*0.8), legend = levels(data$Group), col = code_color, pch = 16)
dev.off()