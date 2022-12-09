library(ggplot2)
library(ggforce)
# source('~/scripts/bin/scale-hue.r')
source('~/scripts/bin/scale-manual.r')

args <- commandArgs(T)
dataFN <- args[1]
outFN <- args[2]
yLabel <- args[3]
yLim <- as.double(strsplit(args[4],split=';')[[1]]) # "yLim_l;yLim_u;breaks" | "0;1;0.5"
order_xValue <- strsplit(args[5],split = ';')[[1]]
order_group <- strsplit(args[6], split = ';')[[1]]
order_col <- strsplit(args[7], split = ';')[[1]]
cols <- strsplit(args[8], split = ';')[[1]] #D9979F;#D95F02;#CE0013;#8DCEBB;#1B9E77;#16557A
jitterPosition <- as.double(args[9])
val_trans <- as.double(args[10])
hLines <- as.double(strsplit(args[11], split = ';')[[1]])
pdf_width <- as.integer(args[12])
pdf_height <- as.integer(args[13])


data <- read.table(dataFN, sep='\t', header=T, stringsAsFactors = F) # format: xValue, yValue, group, col
data <- na.omit(data)
data <- data[data$xValue %in% order_xValue, ]
data$xValue <- factor(data$xValue, levels = order_xValue)
data <- data[data$group %in% order_group, ]
data$group <- factor(data$group, levels = order_group)

# col to numeric
for(i in 1:length(order_col)){
	data[data == order_col[i]] <- i
}
data$col <- as.numeric(data$col)

print(cols)
print(order_col)
print(table(data$col))
pdf(outFN,width=pdf_width,height=pdf_height)
		ggplot(data, aes(xValue, yValue)) +
		# geom_violin() +
		geom_sina(aes(colour = col, shape = group), alpha = val_trans, position = 'identity', seed = 5, scale = F) +
		# geom_text(aes(label = label)) +
		# ylim(c(yLim[1], yLim[2])) +
		scale_colour_gradientn(colours = cols, labels = order_col, breaks = seq(1, length(order_col), 1)) +
		labs(y = yLabel) +
		theme(panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		geom_hline(yintercept = hLines, linetype = 'dashed', color = 'gray')
dev.off()