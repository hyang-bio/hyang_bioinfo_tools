# Usag: 
# Func: quality control of RNA-Seq using hisat
# Date: Mar-01-2019


args <- commandArgs(T)
dirPATH <- args[1] # dirPATH
filePATTEN <- args[2]
outF <- args[3] # outF prefix
cols <- strsplit(args[4], ";")[[1]] # colors for PCA plot


library(gplots)

# combine all fpkm
read.fpkm <- function(dirPATH){
  filename<-dir(dirPATH, pattern = filePATTEN)
  data_first <- read.table(paste(dirPATH,"/",filename[1],sep=""),header=T)
  cln <- c('chr','start','end','strand','t_name','gene_id','t_id','FPKM')
  res <- data_first[,match(cln,colnames(data_first))]
  cln[1] <- '#chrom';cln[8] <- strsplit(filename[1],".ctab")[[1]][1]
  for(i in 2:length(filename)){
    data<-read.table(paste(dirPATH,"/",filename[i],sep=""),header=T)
    res<-cbind(res, data[,'FPKM'])
    cln<-c(cln, strsplit(filename[i],".ctab")[[1]][1])
  }
  colnames(res)<-cln
  return(res)
}


fpkm <- read.fpkm(dirPATH)
fpkm.clean <- log(fpkm[,8:dim(fpkm)[2]] + 1, 2)
write.table(fpkm, file = paste(outF,'.FPKM.merged',sep=''),sep = '\t',quote = F, row.names = F, col.names = T)

# correlation
pdf(paste(outF,'.correlation.heatmap.pdf', sep = ''), width = 7, height = 7, , useDingbats = F)
# par(mar=c(5.1,4.1,0.8,5.1))
cor.val <- cor(fpkm.clean)
cor.lab <- matrix(as.character(round(cor.val,2)),ncol=dim(cor.val)[2])
heatmap.2(cor.val,col=bluered(100),Rowv=T,Colv="Rowv",dendrogram="both",trace="none",cellnote=cor.lab,notecol="black",notecex = 0.5,cexRow=0.5,cexCol=0.5,main="correlation of samples", margins = c(8,8))
dev.off()

pdf(paste(outF,'.correlation.heatmap_primaryOrder.pdf', sep = ''), width = 7, height = 7, , useDingbats = F)
	heatmap.2(cor.val,col=bluered(100),Rowv = F,Colv = F,dendrogram = "none", trace="none",cellnote=cor.lab,notecol="black",notecex = 0.5,cexRow=0.5,cexCol=0.5,main="correlation of samples", margins = c(8,8))
dev.off()

# pca
pdf(paste(outF,'.pca.pdf',sep=''), useDingbats = F)
pc1<-prcomp(t(fpkm.clean))
plot(pc1$x[,1], pc1$x[,2], pch=20, main = "PCA of used samples", xlab = "PC 1", ylab = "PC 2", col = cols)
text(pc1$x[,1], pc1$x[,2],labels=colnames(fpkm.clean),cex=0.8, col = cols)
dev.off()