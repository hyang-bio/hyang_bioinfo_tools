args <- commandArgs(T)

outF <- args[1]
image_width <- as.double(args[2])
image_height <- as.double(args[3])
ratio_height <- as.integer(strsplit(args[4],split=";")[[1]])
zmax <- as.double(args[5])
zmin <- as.double(args[6])
zstep <- as.double(args[7])
fileName <- args[8]
headerLogic <- as.logical(args[9])
cols <- as.integer(strsplit(args[10],split="\n")[[1]]) # 
xlabList <- strsplit(args[11], split = ";")[[1]]
titleList <- args[12]
hLines <- as.double(strsplit(args[13], split = ';')[[1]])
vLines <- as.double(strsplit(args[14], split = ';')[[1]])
format <- args[15]
colors <- args[16] # #DA404E (DNA methylation), #507f9B (H3K9me3)
dataOrder <- args[17]


library('gplots')
library('RColorBrewer')


heatmap_methyl_withRowLabel <- function(matrix1,zmax,zmin,zstep,xlabList, ylab1,Title, hLines, colors){
    matrix1_raw <- matrix1
    matrix1[ matrix1 > zmax ] <- zmax
    matrix1[ matrix1 < zmin ] <- zmin
    matrix1[is.na(matrix1)] <- zmax + zstep # specific color (gray) for NA

    ColorRamp <- colorRampPalette(c('white', colors))(100)
    ColorLevels <- seq(to = zmax, from = zmin, length = 100)
    
   
    # graph 1
    par(mar=c(2,2,6,10)/2)
    image(1:ncol(matrix1), 1:nrow(matrix1), t(matrix1), xaxt = "n", col = c(ColorRamp, 'gray'), xlab = '', ylab = '', axes = FALSE, breaks = c(seq(zmin, zmax, (zmax - zmin)/100), zmax + zstep), zlim = c(zmin, zmax + zstep))
    title(main=Title, cex.main = 0.75)
    axis(side=1,at=seq(from=1,to=ncol(matrix1),1), labels = xlabList,tick=FALSE,las=2, cex.axis = 0.4)
    axis(side=4,at=seq(from=1,to=nrow(matrix1),1),labels=ylab1,tick=FALSE,las=1, cex.axis = 0.3)
    if(hLines[1] != -007){
        for(i in 1:length(hLines)){
            abline(h = hLines[i], lwd = 2, lty = 2)
        }
    }
    if(vLines[1] != -007){
        for(i in 1:length(vLines)){
            abline(v = vLines[i], lwd = 2)
        }
    }

    box()


    # axis
    par(mar=c(4,2,4,10)/2)

    # bottom label
    image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)
    axis(side=1,seq(zmin,zmax,zstep),labels=round(seq(zmin,zmax,zstep),2), cex.axis = 0.8)
    box()
}

heatmap_methyl_withoutRowLabel <- function(matrix1,zmax,zmin,zstep, xlabList, Title, hLines, colors){
    matrix1_raw <- matrix1
    matrix1[ matrix1 > zmax ] <- zmax
    matrix1[ matrix1 < zmin ] <- zmin

    ColorRamp <- colorRampPalette(c('white', colors))(100)
    ColorLevels <- seq(to = zmax, from = zmin, length = 100)
    
   
    # graph 1
    par(mar=c(2,2,3.6,2)/2)
    image(1:ncol(matrix1), 1:nrow(matrix1), t(matrix1), xaxt = "n", col = ColorRamp, xlab = '', ylab = '', axes = FALSE, breaks = seq(zmin, zmax, (zmax - zmin)/100), zlim = c(zmin, zmax))
    title(main=Title, cex.main = 0.75)
    axis(side=1,at = seq(1, ncol(matrix1), 1), labels = xlabList ,las=2, cex.axis = 0.4)
    if(hLines[1] != -007){
        for(i in 1:length(hLines)){
            abline(h = hLines[i], lwd = 2, lty = 2)
        }
    }

    if(vLines[1] != -007){
        for(i in 1:length(vLines)){
            abline(v = vLines[i], lwd = 2)
        }
    }
    box()


    # axis
    par(mar=c(3.7,2,3,2)/2)

    # bottom label
    image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",xaxt="n",yaxt="n",useRaster=T)
    axis(side=1,seq(zmin,zmax,zstep),labels=round(seq(zmin,zmax,zstep),2), cex.axis = 0.4)
    box()
}

if(format == 'pdf'){
    pdf(paste(outF, '.pdf', sep = ''), height = image_height, width = image_width)
}else if(format == 'png'){
    png(paste(outF, '.png', sep = ''), height = image_height, width = image_width, units = 'px', pointsize = 12)
}

layout(matrix(seq(from = 1, to = 2, by = 1),nrow = 2, byrow = F),height = ratio_height)
data1 <- read.table(fileName, sep='\t', header = headerLogic, comment.char = '', check.names = F)
if(dataOrder == "Reverse"){
    data1 <- data1[seq(dim(data1)[1],1), ]
}
if(length(args) == 18){
yCol <- as.integer(args[18])
heatmap_methyl_withRowLabel(as.matrix(data1[,match(cols,seq(1,ncol(data1),1))],ncol = length(cols)),zmax,zmin,zstep, xlabList,as.character(data1[,yCol]),titleList, hLines, colors)
}else if(length(args) == 17){
heatmap_methyl_withoutRowLabel(as.matrix(data1[,match(cols,seq(1,ncol(data1),1))],ncol = length(cols)),zmax,zmin,zstep, xlabList,titleList, hLines, colors)
}
dev.off()