args <- commandArgs(T)
xData <- args[1]
xIndex <- as.integer(args[2])
xLabel <- args[3]
yData <- args[4]
yIndex <- as.integer(args[5])
yLabel <- args[6]
outF <- args[7]
logicL <- args[8]
quanL <- args[9]
xLim <- as.double(strsplit(args[10], split = ';')[[1]]) # 0;1
yLim <- as.double(strsplit(args[11], split = ';')[[1]]) # 0;1
# print(args)


SmoothScatter <- function(xData,xIndex,xLabel,yData,yIndex,yLabel,outF, logicL, quanL){
	if(logicL == 'T' | logicL == 'True' | logicL == 'TRUE'){
		xDataValue <- log2(read.table(xData,sep='\t', header = T)[,xIndex]+1)
		yDataValue <- log2(read.table(yData,sep='\t', header = T)[,yIndex]+1)
	}else if(logicL == 'F' | logicL == 'False' | logicL == 'FALSE'){
		xDataValue <- read.table(xData,sep='\t', header = T)[,xIndex]
		yDataValue <- read.table(yData,sep='\t', header = T)[,yIndex]
	}else{
		print("Need Logic for logicL")
	}

	
	data <- data.frame(xDataValue=xDataValue,yDataValue=yDataValue)
	data <- na.omit(data)

	if(quanL == 'T' | quanL == 'True' | quanL == 'TRUE'){
		quan_value <- quantile(unlist(data), 0.995)
		data[data>quan_value] <- quan_value
	}

	corValue_PE <- cor(data[,1],data[,2], method = 'pearson')
	corValue_SP <- cor(data[,1],data[,2], method = 'spearman')
	corValue_KE <- cor(data[,1],data[,2], method = 'kendall')
	pdf(paste('SmoothScatter of ',outF,'.pdf',sep=''),width=6,height=6)
	par(pty='s')
	smoothScatter(data$xDataValue,data$yDataValue,xlab=xLabel,ylab=yLabel,xlim = xLim, ylim = yLim)
	# add reference line : y = x
	# abline(0,1,col='darkgrey',lty=5,lwd=2)

	# add regression fit line
	fit <- lm(yDataValue~xDataValue,data)
	abline(fit,col='red',lwd=2)
	
	# add correlation value
	text_x <- 0.2*(xLim[2] - xLim[1])
	text_y <- 0.4*(yLim[2] - yLim[1])
	text_y2 <- 0.6*(yLim[2] - yLim[1])
	text_y3 <- 0.8*(yLim[2] - yLim[1])
	text(text_x,text_y,paste('R(PE) = ',round(corValue_PE,2),sep=''))
	text(text_x,text_y2,paste('R(SP) = ',round(corValue_SP,2),sep=''))
	text(text_x,text_y3,paste('R(KE) = ',round(corValue_KE,2),sep=''))
	dev.off()
}
SmoothScatter(xData,xIndex,xLabel,yData,yIndex,yLabel,outF,logicL, quanL)
