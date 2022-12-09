args <- commandArgs(T)
dataFN <- args[1]
outFN <- args[2]
mainLabel <- args[3]
xLabel <- args[4]
yLabel <- args[5]
ylim_l <- as.double(args[6])
ylim_u <- as.double(args[7])
breaks <- as.double(args[8])
order_xValue <- strsplit(args[9],split = ';')[[1]]
order_group <- strsplit(args[10],split=';')[[1]]
pdf_width <- as.double(args[11])
pdf_height <- as.double(args[12])
hLines <- as.double(strsplit(args[13], split = ';')[[1]])

library(ggplot2)
data <- read.table(dataFN,sep='\t',header=T, check.names = F, comment.char = '') # format: xValue, yValue, group
data <- na.omit(data)
data <- data[data$group %in% order_group,]
data <- data[data$xValue %in% order_xValue,]
data$group <- factor(data$group, levels=order_group)
data$xValue <- factor(data$xValue, levels=order_xValue)
pdf(outFN,width=pdf_width,height=pdf_height, useDingbats = F)

p <- ggplot(data,aes(x=xValue,y=yValue,color = group)) +
	geom_boxplot(color = 'black', outlier.colour = NA) +
	geom_jitter(shape = 20, position = position_jitter(0.2)) +
	scale_y_continuous(limits = c(ylim_l,ylim_u), breaks = seq(ylim_l, ylim_u, breaks)) +
	labs(title = mainLabel,x = xLabel,y = yLabel) +
	theme(plot.title = element_text(hjust = 0.5, size = 12)) +
	theme(panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
	theme(axis.text.x = element_text(angle = 90, hjust=1)) +
	geom_hline(yintercept = hLines, linetype = 'dashed', color = 'gray')
	
print(p)
dev.off()