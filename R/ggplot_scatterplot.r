args <- commandArgs(T)
dataFN <- args[1]
mainLabel <- args[2]
xLabel <- args[3]
yLabel <- args[4]
xRange <- as.double(strsplit(args[5], split = ';')[[1]])
yRange <- as.double(strsplit(args[6], split = ';')[[1]])
order_group <- strsplit(args[7], split = ';')[[1]]
Logic_line <- as.logical(strsplit(args[8], split = ';')[[1]][1])
Line_I <- as.double(strsplit(args[8], split = ';')[[1]][2])
Line_S <- as.double(strsplit(args[8], split = ';')[[1]][3])
val_trans <- as.double(args[9])
pdf_width <- as.double(args[10])
pdf_height <- as.double(args[11])
code_color <- args[12]
point_highlight <- args[13]


library(ggplot2)
library(ggExtra)
library(RColorBrewer)
# print(brewer.pal(code_color, n = 8))
data <- read.table(dataFN, sep = '\t', header = T) # format: xValue, yValue, group, shape, ...
data <- na.omit(data)
data <- data[data$group %in% order_group, ]
data$group <- factor(data$group,levels=order_group)

# highlights
highlights <- data.frame(xValue = -007, yValue = -007, group = -007, shape = -007, label = '-007')
if(file.exists(point_highlight)){
	if(file.size(point_highlight) > 0){
		highlights <- read.table(point_highlight, sep = '\t', header = T)
	}
}

# correlation between xValue and yValue
cor_PE <- round(cor(data$xValue, data$yValue), 2)

pdf(paste0('Scatterplot of ', mainLabel, '.pdf'), width = pdf_width,height = pdf_height, useDingbats = F)
if(Logic_line){
	p <- ggplot(data, aes(x = xValue, y = yValue)) +
		geom_point(color = 'gray', alpha = val_trans) + 
		stat_density_2d(aes(fill = ..level..), geom = 'polygon', colour = 'black', size = 0.05) +
		# scale_fill_distiller(palette = 'Purples', directon = 1) +
		scale_fill_gradient(low = '#F9F9FA', high = '#6C3D9A') +
		# geom_density_2d_filled(alpha = val_trans, bins = 5) +
		# geom_density_2d(size = 0.25, colour = 'black', bins = 5) +
		# scale_fill_brewer() +
		# scale_fill_manual(values = colorRampPalette(brewer.pal('Purples', n = 9))(20)) +
		geom_line(intercept = Line_I , slope = Line_S) +
		geom_point(data = highlights, color = '#CE0013', fill = '#CE0013') +
		geom_text(data = highlights, aes(label = label), color = '#CE0013') +
		geom_text(aes(x = xRange[2]*0.2, y = yRange[2]*0.6, label = paste0('r = ', cor_PE))) +
		theme(legend.position = 'left', aspect.ratio = 1) +
		theme(panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
		scale_x_continuous(limits = xRange[1:2], breaks = seq(xRange[1], xRange[2], by = xRange[3]), expand = c(0, 0)) +
		scale_y_continuous(limits = yRange[1:2], breaks = seq(yRange[1], yRange[2], by = yRange[3]), expand = c(0, 0)) +
		labs(title = mainLabel, x = xLabel, y = yLabel) +
		theme(plot.title = element_text(size = 11))
	ggExtra::ggMarginal(p, xparams = list(color = 'black', fill = '#D95F02'), yparams = list(color = 'black', fill = '#1B9E77'), type = 'densigram', alpha = val_trans)
}else{
	p <- ggplot(data, aes(x = xValue, y = yValue)) +
		geom_point(color = 'gray', alpha = val_trans) + 
		stat_density_2d(aes(fill = ..level..), geom = 'polygon', colour = 'black', size = 0.05) +
		scale_fill_gradient(low = '#F9F9FA', high = '#6C3D9A') +
		geom_point(data = highlights, color = '#CE0013', fill = '#CE0013') +
		geom_text(data = highlights, aes(label = label), color = '#CE0013') +
		geom_text(aes(x = xRange[2]*0.2, y = yRange[2]*0.6, label = paste0('r = ', cor_PE))) +
		theme(legend.position = 'left', aspect.ratio = 1) +
		theme(panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
		scale_x_continuous(limits = xRange[1:2], breaks = seq(xRange[1], xRange[2], by = xRange[3])) +
		scale_y_continuous(limits = yRange[1:2], breaks = seq(yRange[1], yRange[2], by = yRange[3])) +
		labs(title = mainLabel, x = xLabel, y = yLabel) +
		theme(plot.title = element_text(size = 11))
	ggExtra::ggMarginal(p, xparams = list(color = 'black', fill = '#D95F02'), yparams = list(color = 'black', fill = '#1B9E77'), type = 'densigram', alpha = val_trans)
}
dev.off()

# ref: https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2.html
# ref: https://biostats.w.uib.no/creating-a-2d-density-plot/