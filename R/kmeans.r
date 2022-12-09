args <- commandArgs(T)
file <- args[1]
headerL <- as.logical(args[2])
colList <- as.double(strsplit(args[3], '\n')[[1]]) # $(seq 4 17)
nCenters <- as.integer(args[4])
algorithmM <- args[5] # "Hartigan-Wong" (not recommend if dots are too close), "Lloyd", "Forgy", "MacQueen"
outFile <- args[6]
wrwnL <- as.logical(args[7])
wclnL <- as.logical(args[8])
# print(colList)

kmeans_ordered <- function(file, headerL, colList, nCenters, outFile, wrwnL, wclnL){
	# read data
	data <- read.table(file, sep = '\t', header = headerL, check.names = F, comment.char = '')
	data <- na.omit(data)

	# kmeans firstly
	set.seed(222)
	ks_1st <- kmeans(data[, colList], centers = nCenters, algorithm = algorithmM, iter.max = 100)

	# kmeans secondly according to specific order
	centerOrder <- ks_1st$centers[order(rowMeans(ks_1st$centers), decreasing = T), ]
	set.seed(222)
	ks_2nd <- kmeans(data[, colList], centers = centerOrder, algorithm = algorithmM, iter.max = 100)

	# write to file
	data_ordered <- data[order(ks_2nd$cluster, decreasing = F), ]
	write.table(data_ordered, file = outFile, sep = '\t', quote = F, row.names = wrwnL, col.names = wclnL)
	write.table(cumsum(rev(ks_2nd$size)), file = paste0(outFile, '.hLines'), sep = '\t', quote = F, row.names = F, col.names = F)

	for(c in 1:nCenters){
		write.table(data[ks_2nd$cluster == c, ], file = paste0('cluster', c, '.txt'), sep = '\t', quote = F, row.names = F, col.names = F)
	}
}

kmeans_ordered(file, headerL, colList, nCenters, outFile, wrwnL, wclnL)