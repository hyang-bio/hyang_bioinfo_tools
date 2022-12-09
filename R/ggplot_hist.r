args <- commandArgs(T)
dataFN <- args[1]
outFN <- args[2]
order_group <- strsplit(args[3],split=';')[[1]]
pdf_width <- as.integer(args[4])
pdf_height <- as.integer(args[5])


library(ggplot2)
data <- read.table(dataFN, sep = '\t', header = T) # format [value group]
data <- na.omit(data)
data <- data[data$group %in% order_group,] 
data$group <- factor(data$group,levels=order_group)

pdf(outFN, width = pdf_width, height = pdf_height)
ggplot(data, aes(x = value, color = group)) +
	geom_histogram(fill = 'white', alpha = 0.5, position = 'identity') +
	# geom_histogram(aes(y = ..density..), binwidth = 1, position = 'dodge') +
	geom_density() + 
	theme(panel.background = element_blank(),panel.border = element_rect(colour = 'black', fill = NA, size = 1))
dev.off()