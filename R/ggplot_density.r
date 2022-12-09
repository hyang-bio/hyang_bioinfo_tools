args <- commandArgs(T)
dataFN <- args[1]
outFN <- args[2]
order_group <- strsplit(args[3],split=';')[[1]]
xlim_l <- as.double(args[4])
xlim_u <- as.double(args[5])
breaks <- as.double(args[6])
pdf_width <- as.integer(args[7])
pdf_height <- as.integer(args[8])


library(ggplot2)
data <- read.table(dataFN, sep = '\t', header = T) # format [value group]
data <- na.omit(data)
data <- data[data$group %in% order_group,] 
data$group <- factor(data$group,levels=order_group)

pdf(outFN, width = pdf_width, height = pdf_height)
ggplot(data, aes(x = value, color = group)) +
	geom_density(alpha = 0.4) + 
	scale_x_continuous(limits = c(xlim_l, xlim_u), breaks = seq(xlim_l, xlim_u, breaks)) +
	theme(panel.background = element_blank(),panel.border = element_rect(colour = 'black', fill = NA, size = 1))
dev.off()