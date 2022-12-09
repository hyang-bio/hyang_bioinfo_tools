args <- commandArgs(T)

library(ggplot2)
library(ggalluvial)

dataFN <- args[1]
outFN <- args[2]
mainLabel <- args[3]
order_group <- strsplit(args[4],split=';')[[1]]
labelS <- strsplit(args[5],split=';')[[1]]
pdf_width <- as.double(args[6])
pdf_height <- as.double(args[7])
cols <- strsplit(args[8], split = ';')[[1]] # "#7570B3;#D95F02;#1B9E77;#CE0013;#16557A"

# dataFN <- 'methylChanges_1kb.groups.freq.txt'
# outFN <- "methylChanges_1kb_Freq_alluvial.pdf"
# mainLabel <- "DNA methylation changes during \nmouse early embryogenesis[1kb]"
# order_group <- strsplit("[0,0.2);[0.2,0.4);[0.4,0.6);[0.6,0.8);[0.8,1.0]",split=";")[[1]]

data <- read.table(dataFN, sep='\t', header=T, check.names = F)
data$From <- factor(data$From,levels = order_group)
data$To <- factor(data$To,levels = order_group)
# all stages together
pdf(outFN, width = pdf_width, height = pdf_height)
ggplot(data,aes(y = Freq, axis1 = From, axis2 = To)) +
  geom_alluvium(aes(fill = From), width = 0, knot.pos = 0, reverse = FALSE) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE, alpha = 1) +
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE) +
  scale_x_continuous(breaks = 1:2, labels = labelS) +
  scale_fill_manual(values = cols) +
  ggtitle(mainLabel)
dev.off()
