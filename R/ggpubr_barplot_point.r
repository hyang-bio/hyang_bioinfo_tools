args <- commandArgs(T)
dataN <- args[1]
outFN <- args[2]
yLabel <- gsub('__n', '\n', args[3])
yRange <- as.double(strsplit(args[4], split = ';')[[1]])
xValue_order <- strsplit(args[5], split = ';')[[1]]
group_order <- strsplit(args[6], split = ';')[[1]]
pdf_width <- as.double(args[7])
pdf_height <- as.double(args[8])
hLines <- as.double(strsplit(args[9], split = ';')[[1]])
percentageL <- args[10]


library(ggpubr)
value <- read.table(dataN, sep='\t', header = T, check.names = F, comment.char = '') # file format: xValue, yValue, group
value <- value[which(value$xValue %in% xValue_order & value$group %in% group_order),]
value$group <- factor(value$group, levels = group_order)
value$xValue <- factor(value$xValue, levels = xValue_order)
if(percentageL == 'T' | percentageL == 'True' | percentageL == 'TRUE'){
	value$yValue <- value$yValue*100
}

pdf(outFN, width = pdf_width, height = pdf_height, useDingbats = F)
	ggbarplot(value,
		x = 'xValue',
		y = 'yValue',
		xlab = '',
		ylab = yLabel,
		add = c('mean', 'point'),
		fill = 'group',
		palette = 'npg',
		position = position_dodge(0.8)) +
	scale_y_continuous(limits = yRange[1:2], breaks = seq(yRange[1], yRange[2], by = yRange[3])) +
	geom_hline(yintercept = hLines, linetype = 'dashed', color = 'gray') + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()