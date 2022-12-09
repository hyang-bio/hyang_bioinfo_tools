args <- commandArgs(T)
dataN <- args[1]
mainLabel <- args[2]
group_order <- strsplit(args[3], split = ';')[[1]]
cols <- strsplit(args[4], split = ';')[[1]]
percentageL <- args[5]
pdf_width <- as.double(args[6])
pdf_height <- as.double(args[7])


library(ggplot2)
library(RColorBrewer)
library(scales)

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
)

df <- read.table(dataN, sep = '\t', header = T, check.names = F) # column: value, group
df <- df[match(group_order, df$group),]
df$group <- factor(df$group, levels = group_order)
df <- na.omit(df)

if(percentageL == 'T' | percentageL == 'True' | percentageL == 'TRUE'){
	df$value <- df$value/sum(df$value)*100
}

pdf(paste0('Pie plot of ',mainLabel,'.pdf'), width = pdf_width, height = pdf_height)
ggplot(df, aes(x = "", y = value, fill = group)) +
	geom_bar(width = 1, stat = 'identity') +
	coord_polar('y', start = 0) +
	blank_theme +
	theme(axis.text.x = element_blank()) +
	geom_text(aes(y = (100- value/2 - c(0, cumsum(value)[-length(value)])), label = round(value, 2))) +
	# scale_fill_manual(values = colorRampPalette(brewer.pal('Dark2', n = 8))(nlevels(df$group))) +
	scale_fill_manual(values = cols) +
	labs(title=mainLabel)
dev.off()