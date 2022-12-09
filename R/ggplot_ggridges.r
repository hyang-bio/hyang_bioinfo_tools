args <- commandArgs(T)
dataN <- args[1]
outFN <- args[2]
label <- args[3]
group_order <- strsplit(args[4], split = ';')[[1]]
xlim_l <- as.double(args[5])
xlim_u <- as.double(args[6])
breaks <- as.double(args[7])
pdf_width <- as.double(args[8])
pdf_height <- as.double(args[9])

library(ggplot2)
library(ggridges)
library(RColorBrewer)

data <- read.table(dataN, sep = '\t', header = T) # format: value, group
data <- data[data$group %in% group_order, ]
data$group <- factor(data$group, levels = group_order)

pdf(outFN,width=pdf_width,height=pdf_height, useDingbats = F)
	ggplot(data, aes(x = value, y = group)) +
	geom_density_ridges(aes(fill = group), scale = 5, alpha = 0.7) +
	scale_x_continuous(limits = c(xlim_l, xlim_u), breaks = seq(xlim_l, xlim_u, breaks)) +
	labs(x = label) +
	theme(panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
	scale_fill_manual(values = colorRampPalette(brewer.pal('Dark2', n = 8))(nlevels(data$group)))
dev.off()