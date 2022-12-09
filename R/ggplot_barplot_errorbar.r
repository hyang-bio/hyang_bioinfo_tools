args <- commandArgs(T)
dataN <- args[1]
outFN <- args[2]
mainLabel <- args[3]
xLabel <- args[4]
yLabel <- gsub('__n', '\n', args[5])
yRange <- as.double(strsplit(args[6], split = ';')[[1]])
xValue_order <- strsplit(args[7],split=';')[[1]]
group_order <- strsplit(args[8],split=';')[[1]]
position_type <- args[9] # dodge(并列), identity(上下)
pdf_width <- as.double(args[10])
pdf_height <- as.double(args[11])
hLines <- as.double(strsplit(args[12], split = ';')[[1]])
percentageL <- args[13]

library(ggplot2)
library(RColorBrewer)
value <- read.table(dataN, sep='\t', header = T, check.names = F, comment.char = '') # file format: xValue, yValue, sd, group
value <- value[which(value$xValue %in% xValue_order & value$group %in% group_order),]
value$group <- factor(value$group,levels=group_order)
value$xValue <- factor(value$xValue,levels=xValue_order)
if(percentageL == 'T' | percentageL == 'True' | percentageL == 'TRUE'){
	value$yValue <- value$yValue*100
}

# error bar
value$sd_up <- value$yValue + value$sd
value$sd_low <- value$yValue - value$sd

pdf(outFN, width=pdf_width,height=pdf_height)
ggplot(value,aes(x=xValue,y=yValue,fill=group)) +
	geom_errorbar(aes(ymin = sd_low, ymax = sd_up), position = position_type) +
	scale_fill_manual(values = colorRampPalette(brewer.pal('Dark2', n = 8))(nlevels(value$group))) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	scale_y_continuous(limits = yRange[1:2], breaks = seq(yRange[1], yRange[2], by = yRange[3])) +
	labs(title= mainLabel,x = xLabel ,y = yLabel) +
	theme(panel.background = element_blank(),panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
	geom_hline(yintercept = hLines, linetype = 'dashed', color = 'gray')
   
dev.off()