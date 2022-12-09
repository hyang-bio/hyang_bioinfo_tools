args <- commandArgs(T)

dataFN <- args[1]
outFN <- args[2]
mainLabel <- args[3]
order_group <- strsplit(args[4],split=';')[[1]]

library(ggplot2)
library(plotROC)
data <- read.table(dataFN, sep = '\t', header =T, check.names = F, comment.char = '') # column necessary: D(Decision: 0/1) D.str(Decision: Fail/Success) M(Measure: value) group
data$group <- factor(data$group, levels = order_group)

pdf(outFN,width = 6, height =6)
par(pty = "s")
basicplot <- ggplot(data, aes(d = D, m = M, color = group)) + 
geom_roc(n.cuts = 0) +
scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
labs(title=mainLabel,x='False Positive Fraction (1-Specificity)',y='True Positive Fraction (Sensitivity)') +
theme(plot.title = element_text(hjust = 0.5, size = 10)) +
theme(panel.background = element_blank(),panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
theme(aspect.ratio = 1)
direct_label(basicplot, labels = paste(rep("AUC(",2),order_group,rep(")=",2),round(calc_auc(basicplot)$AUC,3),sep=''),label.angle = 0, nudge_y = -.1, size = 3)
res <- c(dataFN, round(calc_auc(basicplot)$AUC,3))
write.table(t(res), file = paste(outFN, '.AUCs', sep = ''), quote = F, row.names = F, col.names = F)
dev.off()