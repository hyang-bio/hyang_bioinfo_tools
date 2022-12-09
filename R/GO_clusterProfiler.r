args <- commandArgs(T)
genesFile <- args[1] 
geneType <- args[2] # REFSEQ
species <- args[3] # Mm
outFilePre <- args[4]
pdf_width <- as.double(args[5])
pdf_height <- as.double(args[6])

library(ggplot2)
library(clusterProfiler)

GO_clusterProfiler <- function(genesFile, geneType, species, outFilePre, pdf_width, pdf_height){
	# step1. convert to ENTREZID
	genes <- as.character(read.table(genesFile, sep = '\t', header = F)$V1)
	eg <- bitr(geneID = genes, fromType = geneType, toType = 'ENTREZID', OrgDb = paste0('org.', species, '.eg.db'))

	# step2. enrich GO
	go <- enrichGO(eg$ENTREZID, 
		OrgDb = paste0('org.', species, '.eg.db'), 
		ont = 'BP', 
		pAdjustMethod = 'BH', 
		pvalueCutoff = 0.01, 
		qvalueCutoff = 0.05)

	# step3. plot
	dotplot(go)
	ggsave(paste0(outFilePre, '.pdf'), width = pdf_width, height = pdf_height)
}

GO_clusterProfiler(genesFile, geneType, species, outFilePre, pdf_width, pdf_height)